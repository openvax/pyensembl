"""Shared canonical Ensembl DNA and conservative orphan pruning.

The versioned assembly accession and file metadata identify an upstream
artifact. Ensembl CHECKSUMS are Unix checksums, not cryptographic digests.
Custom sources are deliberately excluded from this namespace.
"""

from contextlib import contextmanager
import hashlib
import json
import logging
from pathlib import Path
import re
import shutil
from urllib.error import URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

from filelock import FileLock

from .download_cache import DownloadCache
from .genome_fasta import GenomeFasta, _fingerprint, _read_json, _write_json

logger = logging.getLogger(__name__)
_KEY = re.compile(r"^[0-9a-f]{64}$")
_SOURCE_PATH = re.compile(r"^/pub/release-(\d+)/(?:[^/]+/)?fasta/([^/]+)/dna/([^/]+)$")
_PROVIDERS = {"ftp.ensembl.org", "ftp.ensemblgenomes.ebi.ac.uk"}


def dna_cache_root():
    """Shared DNA root under the same global cache as annotation data."""
    return Path(DownloadCache(None, None).cache_directory_path) / "dna_cache"


def is_canonical_source(source):
    url = urlsplit(source)
    return (
        url.scheme == "https"
        and url.netloc in _PROVIDERS
        and _SOURCE_PATH.fullmatch(url.path) is not None
    )


def _identity_key(identity):
    return hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()


def _remote_identity(source):
    """Fetch small metadata; incomplete identity falls back to a single URL."""
    fallback = {"source": source}
    if not is_canonical_source(source):
        return fallback
    url = urlsplit(source)
    release, species, filename = _SOURCE_PATH.fullmatch(url.path).groups()
    directory = source.rsplit("/", 1)[0]
    try:
        with urlopen(directory + "/README", timeout=60) as response:
            readme = response.read(1024 * 1024).decode("utf-8")
        assembly = re.search(r"GenBank Assembly ID\s+(GCA_\d+\.\d+)", readme)
        if assembly is None:
            raise ValueError("No versioned assembly accession")
        with urlopen(directory + "/CHECKSUMS", timeout=60) as response:
            checksums = response.read(4 * 1024 * 1024).decode("utf-8")
        entries = [line.split() for line in checksums.splitlines()]
        checksum = next(
            (
                fields[:2]
                for fields in entries
                if len(fields) == 3
                and fields[2] == filename
                and all(value.isdigit() for value in fields[:2])
            ),
            None,
        )
        if checksum is None:
            raise ValueError("No checksum for this FASTA")
        with urlopen(Request(source, method="HEAD"), timeout=60) as response:
            length = int(response.headers["Content-Length"])
        if assembly is not None and checksum is not None and length > 0:
            return {
                "provider": url.netloc,
                "species": species,
                "assembly": assembly.group(1),
                "filename": filename.replace(".%s.dna" % release, ".dna", 1),
                "checksum": [int(value) for value in checksum],
                "compressed_size": length,
            }
    except (OSError, URLError, ValueError, TypeError, KeyError, UnicodeError):
        pass
    logger.warning(
        "Cannot establish shared DNA identity for %s; using a release-specific object",
        source,
    )
    return fallback


@contextmanager
def dna_cache_lock(root=None):
    root = Path(root) if root is not None else dna_cache_root()
    root.mkdir(parents=True, exist_ok=True)
    with FileLock(str(root / ".lock")):
        yield root


class SharedGenomeFasta(GenomeFasta):
    """Canonical DNA with a per-release reference and a shared immutable key."""

    def __init__(self, source, cache_directory):
        if not is_canonical_source(source):
            raise ValueError("Only official Ensembl DNA can use the shared cache")
        super().__init__(source, cache_directory)
        # Native Ensembl annotations live at cache_root/reference/ensemblN.
        self.root = Path(cache_directory).parent.parent / "dna_cache"
        source_key = hashlib.sha256(self.source.encode()).hexdigest()
        self.reference_path = (
            Path(cache_directory) / "genome_fasta_refs" / (source_key + ".json")
        )
        self.key = None
        self.identity = None
        self._load_reference()

    def _load_reference(self):
        state = _read_json(self.reference_path)
        if (
            isinstance(state, dict)
            and state.get("source") == self.source
            and "shared_key" in state
        ):
            key = state["shared_key"]
            identity = state.get("identity")
            if (
                not isinstance(key, str)
                or not _KEY.fullmatch(key)
                or not isinstance(identity, dict)
                or _identity_key(identity) != key
            ):
                raise ValueError(
                    "Invalid shared genome FASTA reference: %s" % self.reference_path
                )
            self._select(key, identity)

    def _select(self, key, identity):
        if self.key != key:
            self.close()
        self.key = key
        self.identity = identity
        self._use_directory(self.root / "objects" / key)

    def _remember_unlocked(self):
        state = {
            "source": self.source,
            "shared_key": self.key,
            "identity": self.identity,
        }
        _write_json(self.reference_path, state)
        _write_json(self.manifest_path, state)

    def prepare(self, download=False, overwrite=False):
        # A read of an already registered object does not need a writer lock.
        # Its release manifest protects it from pruning.
        if not download:
            return super().prepare()
        with dna_cache_lock(self.root):
            self._load_reference()
            if self.key is None or overwrite:
                identity = _remote_identity(self.source)
                self._select(_identity_key(identity), identity)
            self.directory.mkdir(parents=True, exist_ok=True)
            _write_json(self.directory / "object.json", {"identity": self.identity})
            path = super().prepare(download=True, overwrite=overwrite)
            self._remember_unlocked()
            return path

    def remember(self):
        with dna_cache_lock(self.root):
            if self.key is not None and self.installed_path is not None:
                self._remember_unlocked()

    def _materialize(self, stream):
        super()._materialize(stream, expected_size=self.identity.get("compressed_size"))

    def open(self, overwrite=False):
        if self._reader is not None and not overwrite:
            if self._reader_fingerprint == _fingerprint(self.prepare()):
                return self._reader
        with dna_cache_lock(self.root):
            reader = super().open(overwrite=overwrite)
            self._remember_unlocked()
            return reader


def _reference_manifests(root):
    # Path.glob can silently skip unreadable directories. Pruning must fail
    # closed if any annotation directory cannot be inspected.
    try:
        for reference in root.parent.iterdir():
            if reference == root or not reference.is_dir():
                continue
            for annotation in reference.iterdir():
                if not annotation.is_dir():
                    continue
                manifest = annotation / "genome_fasta.json"
                try:
                    manifest.lstat()
                except FileNotFoundError:
                    pass
                else:
                    yield manifest
                try:
                    for path in (annotation / "genome_fasta_refs").iterdir():
                        if path.name.endswith(".json"):
                            yield path
                except FileNotFoundError:
                    pass
    except OSError as error:
        raise ValueError(
            "Cannot safely prune: unable to inspect release references"
        ) from error


def prune_genome_fastas(dry_run=False, cache_root=None):
    """Remove owned shared objects with no release references.

    Return (path, bytes) pairs. Malformed manifests abort the entire operation;
    dry_run performs the same reference checks without removing anything.
    Attached local files, private downloads, and symlinks are never candidates.
    """
    root = Path(cache_root) if cache_root is not None else dna_cache_root()
    if not root.exists():
        return []
    with dna_cache_lock(root):
        referenced = set()
        for path in _reference_manifests(root):
            state = _read_json(path)
            if not isinstance(state, dict) or not isinstance(state.get("source"), str):
                raise ValueError(
                    "Cannot safely prune: invalid genome FASTA manifest %s" % path
                )
            if "shared_key" in state:
                key = state["shared_key"]
                identity = state.get("identity")
                if (
                    not isinstance(key, str)
                    or not _KEY.fullmatch(key)
                    or not isinstance(identity, dict)
                    or _identity_key(identity) != key
                ):
                    raise ValueError(
                        "Cannot safely prune: invalid shared reference %s" % path
                    )
                referenced.add(key)
        candidates = []
        objects = root / "objects"
        if objects.is_symlink():
            raise ValueError("Refusing to prune a symlinked DNA objects directory")
        if not objects.exists():
            return []
        for directory in sorted(objects.iterdir()):
            if (
                directory.is_symlink()
                or not directory.is_dir()
                or not _KEY.fullmatch(directory.name)
            ):
                continue
            if directory.name in referenced:
                continue
            state = _read_json(directory / "object.json")
            if not isinstance(state, dict) or not isinstance(
                state.get("identity"), dict
            ):
                continue
            if _identity_key(state["identity"]) != directory.name:
                continue
            # Never traverse symlinked entries or report external file sizes.
            size = sum(
                path.lstat().st_size
                for path in directory.iterdir()
                if not path.is_dir()
            )
            candidates.append((str(directory), size))
        if not dry_run:
            for path, _ in candidates:
                shutil.rmtree(path)
        return candidates
