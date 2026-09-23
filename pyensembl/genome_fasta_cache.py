"""Shared canonical Ensembl DNA and conservative orphan pruning.

The versioned assembly accession and file metadata identify an upstream
artifact. Ensembl CHECKSUMS are Unix checksums, not cryptographic digests.
Custom sources are deliberately excluded from this namespace.
"""

from contextlib import contextmanager
import hashlib
import json
import logging
import os
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
_COMPONENT = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")
_ACCESSION = re.compile(r"GCA_\d+\.\d+")
_DNA_FILENAME = re.compile(
    r"(.+)\.(dna|dna_sm|dna_rm)\.(toplevel|primary_assembly)\.fa\.gz"
)


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


def _object_directory(root, identity):
    """Expose biological dimensions in the path; keep file revisions distinct."""
    if "source" in identity:
        source = identity["source"]
        if not isinstance(source, str) or not is_canonical_source(source):
            raise ValueError("Invalid source in shared genome FASTA identity")
        url = urlsplit(source)
        release, species, filename = _SOURCE_PATH.fullmatch(url.path).groups()
        filename = filename.replace(".%s.dna" % release, ".dna", 1)
        provider, accession = url.netloc, "unverified"
    else:
        provider = identity.get("provider")
        species = identity.get("species")
        filename = identity.get("filename")
        accession = identity.get("assembly")
        if not isinstance(accession, str) or not _ACCESSION.fullmatch(accession):
            raise ValueError("Invalid assembly in shared genome FASTA identity")
    if not isinstance(provider, str) or provider not in _PROVIDERS:
        raise ValueError("Invalid provider in shared genome FASTA identity")
    if not isinstance(species, str) or not _COMPONENT.fullmatch(species):
        raise ValueError("Invalid species in shared genome FASTA identity")
    prefix = species.capitalize() + "."
    if not isinstance(filename, str) or not filename.startswith(prefix):
        raise ValueError("Invalid filename in shared genome FASTA identity")
    match = _DNA_FILENAME.fullmatch(filename[len(prefix) :])
    if match is None:
        raise ValueError("Invalid DNA format in shared genome FASTA identity")
    reference, sequence_type, coverage = match.groups()
    if not _COMPONENT.fullmatch(reference):
        raise ValueError("Invalid reference in shared genome FASTA identity")
    masking = {"dna": "unmasked", "dna_sm": "softmasked", "dna_rm": "hardmasked"}[
        sequence_type
    ]
    return (
        Path(root)
        / species
        / provider
        / (reference + "-" + accession)
        / coverage
        / masking
        / "fasta"
        / _identity_key(identity)[:16]
    )


def _check_object_identity(root, directory, identity):
    # Metadata-derived paths may never traverse links out of the owned tree.
    for path in (directory, *directory.parents):
        if path.is_symlink():
            raise ValueError("Symlink in shared genome FASTA path: %s" % path)
        if path == root:
            break
    if not directory.exists():
        return
    if not directory.is_dir() or (directory / "object.json").is_symlink():
        raise ValueError("Invalid shared genome FASTA object: %s" % directory)
    state = _read_json(directory / "object.json")
    if state is None and not any(directory.iterdir()):
        return  # An interrupted mkdir has not published any files yet.
    if not isinstance(state, dict) or state.get("identity") != identity:
        # Keep the full identity authoritative, even if short leaf keys collide.
        raise ValueError("Conflicting shared genome FASTA identity: %s" % directory)


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
        directory = _object_directory(self.root, identity)
        _check_object_identity(self.root, directory, identity)
        if self.key != key:
            self.close()
        self.key = key
        self.identity = identity
        self._use_directory(directory)

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


def _owned_objects(root):
    """Find semantic object leaves without following symlinked directories."""
    if root.is_symlink():
        raise ValueError("Refusing to prune a symlinked DNA cache")

    def inaccessible(error):
        raise ValueError(
            "Cannot safely prune: unable to inspect DNA objects"
        ) from error

    for current, directories, files in os.walk(
        root, followlinks=False, onerror=inaccessible
    ):
        directory = Path(current)
        directories[:] = [
            name for name in directories if not (directory / name).is_symlink()
        ]
        if "object.json" not in files or (directory / "object.json").is_symlink():
            continue
        state = _read_json(directory / "object.json")
        if not isinstance(state, dict) or not isinstance(state.get("identity"), dict):
            continue
        identity = state["identity"]
        try:
            expected = _object_directory(root, identity)
        except ValueError:
            continue
        if directory == expected:
            yield directory, _identity_key(identity)


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
                try:
                    _object_directory(root, identity)
                except ValueError as error:
                    raise ValueError(
                        "Cannot safely prune: invalid shared reference %s" % path
                    ) from error
                referenced.add(key)
        candidates = []
        for directory, key in sorted(_owned_objects(root)):
            if key in referenced:
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
