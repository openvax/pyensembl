"""Disk-backed reference DNA, with indexes owned by PyEnsembl."""

from contextlib import contextmanager
import gzip
import hashlib
import json
import logging
import os
from pathlib import Path
import shutil
from urllib.parse import urlsplit
from uuid import uuid4

from datacache import fetch_file

from .normalization import normalize_chromosome

logger = logging.getLogger(__name__)
_CHUNK_SIZE = 1024 * 1024


class MissingGenomeFastaError(ValueError):
    """Reference DNA was not configured or has not been downloaded."""


def _read_json(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def _staging_path(directory, name):
    # Dot files in cache directories are always staging files.
    return Path(directory) / (".%s.%s.tmp" % (name, uuid4().hex))


def _remove(path):
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass


def _publish(staged, path):
    """Durably replace path with a complete staged file."""
    with open(staged, "r+b") as handle:
        os.fsync(handle.fileno())
    os.replace(staged, path)
    try:
        directory = os.open(os.path.dirname(path), os.O_RDONLY)
    except OSError:
        return  # Directories cannot be opened for fsync on some platforms.
    try:
        os.fsync(directory)
    except OSError:
        pass
    finally:
        os.close(directory)


@contextmanager
def _atomic_output(path, mode="wb"):
    """Write a file that appears complete or not at all.

    The staging file is created like any new file (0o666 less umask), so a
    group-shared cache stays readable by other users.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    staged = _staging_path(path.parent, path.name)
    try:
        descriptor = os.open(staged, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666)
        with open(descriptor, mode) as handle:
            yield handle
        _publish(staged, path)
    finally:
        _remove(staged)


def _write_json(path, value):
    with _atomic_output(path, "w") as handle:
        json.dump(value, handle)


def _fingerprint(path, ctime=True):
    stat = os.stat(path)
    if not ctime:
        return [stat.st_size, stat.st_mtime_ns, stat.st_ino]
    return [stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, stat.st_ino]


class GenomeFasta:
    """A single combined FASTA; never modifies an attached local source.

    Compressed input is streamed into an uncompressed, atomically published
    cache file. pyfaidx then reads only the requested interval from disk.
    """

    def __init__(self, source, cache_directory, install_string_function=None):
        self.source = os.fspath(source)
        if not self.source:
            raise ValueError("Genome FASTA source must not be empty")
        self.remote = "://" in self.source
        if not self.remote:
            self.source = os.path.abspath(self.source)
        self.install_string_function = install_string_function
        source_key = hashlib.sha256(self.source.encode()).hexdigest()
        self._use_directory(Path(cache_directory) / "genome_fasta" / source_key)
        self.manifest_path = Path(cache_directory) / "genome_fasta.json"
        self._forget_reader()

    def _use_directory(self, directory):
        self.directory = directory
        self.materialized_path = self.directory / "sequence.fa"
        self.index_path = self.directory / "sequence.fa.fai"
        self.source_state_path = self.directory / "source.json"
        self.index_state_path = self.directory / "index.json"

    def _file_fingerprint(self, path):
        # Cache-owned files change only by atomic replacement (a new inode).
        # Ignore their ctime, which chmod/chown of a shared cache also update;
        # keep it for user files, which may be rewritten in place.
        return _fingerprint(path, ctime=path == self.source)

    def _file_fingerprints(self, paths):
        try:
            return [self._file_fingerprint(path) for path in paths]
        except OSError:
            return None

    def _local_is_compressed(self):
        with open(self.source, "rb") as handle:
            return handle.read(2) == b"\x1f\x8b"

    @property
    def expected_path(self):
        if self.remote or (os.path.isfile(self.source) and self._local_is_compressed()):
            return str(self.materialized_path)
        return self.source

    @property
    def installed_path(self):
        path = self.expected_path
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            return None
        if not self.remote and path != self.source:
            if _read_json(self.source_state_path) != _fingerprint(self.source):
                return None
        return path

    def _not_installed(self):
        message = "Reference DNA is not installed: %s. Call download_genome_fasta()" % (
            self.source
        )
        if self.install_string_function is not None:
            message += " or run: %s" % self.install_string_function()
        return MissingGenomeFastaError(message)

    def _write_uncompressed(self, source_path, destination):
        """Publish plain FASTA from a plain or gzip file, or nothing on failure."""
        with open(source_path, "rb") as raw:
            compressed = raw.peek(2)[:2] == b"\x1f\x8b"
            reader = gzip.GzipFile(fileobj=raw) if compressed else raw
            with reader:
                first = reader.read(_CHUNK_SIZE)
                if not first.startswith(b">"):
                    raise ValueError(
                        "Genome FASTA must start with a FASTA header: %s" % self.source
                    )
                with _atomic_output(destination) as output:
                    output.write(first)
                    shutil.copyfileobj(reader, output, length=_CHUNK_SIZE)

    def _download(self, expected_size=None):
        """Fetch with datacache retries, then publish the decompressed FASTA.

        expected_size describes the downloaded (possibly compressed) bytes.
        """
        self.directory.mkdir(parents=True, exist_ok=True)
        # Keep the URL's suffix so datacache stores the bytes unchanged.
        name = os.path.basename(urlsplit(self.source).path) or "download"
        raw = self.directory / (".download-%s-%s" % (uuid4().hex, name))
        try:
            logger.info("Downloading genome FASTA %s", self.source)
            fetch_file(
                self.source,
                destination=raw,
                force=True,
                timeout=3600,
                expected_size=expected_size,
            )
            self._write_uncompressed(raw, self.materialized_path)
        finally:
            _remove(raw)

    def prepare(self, download=False, overwrite=False):
        """Resolve local data; only an explicit download may access the network."""
        if self.remote:
            if overwrite or self.installed_path is None:
                if not download:
                    raise self._not_installed()
                self._download()
        else:
            if not os.path.isfile(self.source):
                raise MissingGenomeFastaError(
                    "Missing local genome FASTA: %s" % self.source
                )
            if self._local_is_compressed() and (
                overwrite or self.installed_path is None
            ):
                logger.info("Decompressing genome FASTA %s", self.source)
                fingerprint = _fingerprint(self.source)
                self._write_uncompressed(self.source, self.materialized_path)
                if fingerprint != _fingerprint(self.source):
                    raise ValueError(
                        "Genome FASTA changed during decompression: %s" % self.source
                    )
                _write_json(self.source_state_path, fingerprint)
        path = self.installed_path
        if path is None:
            raise self._not_installed() if self.remote else MissingGenomeFastaError(
                "Empty or missing genome FASTA: %s" % self.source
            )
        return path

    def remember(self):
        """Record the attachment for CLI inspection, without changing defaults."""
        state = {"source": self.source}
        if _read_json(self.manifest_path) != state:
            _write_json(self.manifest_path, state)

    def _index_is_current(self, fingerprint):
        return self.index_path.exists() and _read_json(self.index_state_path) == fingerprint

    def _ensure_index(self, path, fingerprint, overwrite=False):
        """Publish a complete index, keeping the old one if parsing fails."""
        from pyfaidx import Fasta

        if not overwrite and self._index_is_current(fingerprint):
            return
        logger.info("Indexing genome FASTA %s", path)
        self.directory.mkdir(parents=True, exist_ok=True)
        staged = _staging_path(self.directory, self.index_path.name)
        try:
            with Fasta(path, indexname=str(staged), strict_bounds=True) as fasta:
                if not len(fasta.keys()):
                    raise ValueError(
                        "Genome FASTA contains no sequences: %s" % self.source
                    )
            if fingerprint != self._file_fingerprint(path):
                raise ValueError("Genome FASTA changed during indexing: %s" % self.source)
            _publish(staged, self.index_path)
            _write_json(self.index_state_path, fingerprint)
        finally:
            _remove(staged)

    def open(self, overwrite=False):
        """Return a pyfaidx reader, building or refreshing the cached index.

        A reader stays in use until its files change; checking that costs a
        stat per file, with no JSON reads, locks, or writes.
        """
        from pyfaidx import Fasta

        if (
            self._reader is not None
            and not overwrite
            and self._file_fingerprints(self._reader_paths) == self._reader_fingerprints
        ):
            return self._reader
        # Callers may still hold the previous reader. Atomically replaced files
        # keep its inode valid, and pyfaidx closes it once unreferenced.
        self._forget_reader()
        path = self.prepare()
        paths = [path] if path == self.source or self.remote else [path, self.source]
        fingerprints = self._file_fingerprints(paths)
        if fingerprints is None:
            raise MissingGenomeFastaError("Genome FASTA disappeared: %s" % path)
        self._ensure_index(path, fingerprints[0], overwrite=overwrite)
        self._reader = Fasta(
            path,
            indexname=str(self.index_path),
            strict_bounds=True,
            build_index=False,
            rebuild=False,
        )
        self._reader_paths = paths
        self._reader_fingerprints = fingerprints
        return self._reader

    def record(self, contig):
        """The pyfaidx record for a contig, or None if absent.

        Accepts the FASTA's own names and pyensembl's normalized contig names:
        annotations store e.g. ``Pt`` as ``PT`` and ``Mito`` as ``MITO``.
        """
        fasta = self.open()
        name = str(contig)
        if name in fasta:
            return fasta[name]
        if self._names_by_normalized is None:
            self._names_by_normalized = {}
            for record in fasta.keys():
                try:
                    key = normalize_chromosome(record)
                except (TypeError, ValueError):
                    continue
                self._names_by_normalized.setdefault(key, []).append(record)
        try:
            matches = self._names_by_normalized.get(normalize_chromosome(contig), [])
        except (TypeError, ValueError):
            return None
        if len(matches) > 1:
            raise ValueError(
                "Contig %r matches several FASTA records: %s" % (contig, ", ".join(matches))
            )
        return fasta[matches[0]] if matches else None

    def _forget_reader(self):
        self._reader = None
        self._reader_paths = ()
        self._reader_fingerprints = None
        self._names_by_normalized = None

    def clear_cache(self):
        """Stop reusing the current reader without closing it for other holders."""
        self._forget_reader()

    def close(self):
        if self._reader is not None:
            self._reader.close()
        self._forget_reader()

    @classmethod
    def installed_source(cls, cache_directory):
        """The DNA recorded for a cache directory, or None.

        Raises ValueError for an invalid shared reference.
        """
        state = _read_json(Path(cache_directory) / "genome_fasta.json")
        if isinstance(state, dict) and isinstance(state.get("source"), str):
            if "shared_key" in state:
                from .genome_fasta_cache import SharedGenomeFasta

                return SharedGenomeFasta(state["source"], cache_directory)
            return cls(state["source"], cache_directory)
        return None

    def status(self, check=False):
        path = self.installed_path
        origin = "downloaded" if self.remote else "local"
        if path is None:
            return "%s, missing" % origin
        indexed = self._index_is_current(self._file_fingerprint(path))
        if check and indexed:
            from pyfaidx import Fasta

            try:
                with Fasta(
                    path,
                    indexname=str(self.index_path),
                    build_index=False,
                    rebuild=False,
                    strict_bounds=True,
                ) as fasta:
                    if not len(fasta.keys()):
                        raise ValueError("empty index")
                    for record in fasta:
                        if len(record) and len(record[-1:].seq) != 1:
                            raise ValueError("truncated FASTA")
            except (OSError, ValueError, IndexError, KeyError):
                return "%s, invalid index: %s" % (origin, path)
        return "%s, %s: %s" % (origin, "indexed" if indexed else "needs index", path)
