"""Disk-backed reference DNA, with indexes owned by PyEnsembl."""

import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import tempfile
from urllib.request import urlopen


class MissingGenomeFastaError(ValueError):
    """Reference DNA was not configured or has not been downloaded."""


class _CountingReader(io.RawIOBase):
    def __init__(self, source):
        self.source = source
        self.bytes_read = 0

    def readable(self):
        return True

    def readinto(self, buffer):
        data = self.source.read(len(buffer))
        buffer[: len(data)] = data
        self.bytes_read += len(data)
        return len(data)


def _read_json(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def _write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
        temporary = handle.name
        try:
            json.dump(value, handle)
            handle.close()
            os.replace(temporary, path)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)


def _fingerprint(path):
    stat = os.stat(path)
    return [stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, stat.st_ino]


class GenomeFasta:
    """A single combined FASTA; never modifies an attached local source.

    Compressed input is streamed into an uncompressed, atomically published
    cache file. pyfaidx then reads only the requested interval from disk.
    """

    def __init__(self, source, cache_directory):
        self.source = os.fspath(source)
        if not self.source:
            raise ValueError("Genome FASTA source must not be empty")
        self.remote = "://" in self.source
        if not self.remote:
            self.source = os.path.abspath(self.source)
        source_key = hashlib.sha256(self.source.encode()).hexdigest()
        self._use_directory(Path(cache_directory) / "genome_fasta" / source_key)
        self.manifest_path = Path(cache_directory) / "genome_fasta.json"
        self._reader = None
        self._reader_fingerprint = None

    def _use_directory(self, directory):
        self.directory = directory
        self.materialized_path = self.directory / "sequence.fa"
        self.index_path = self.directory / "sequence.fa.fai"
        self.source_state_path = self.directory / "source.json"
        self.index_state_path = self.directory / "index.json"

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

    def _materialize(self, stream, expected_size=None):
        self.directory.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=self.directory, delete=False) as output:
            temporary = output.name
            try:
                counted = _CountingReader(stream)
                with io.BufferedReader(counted) as buffered:
                    if buffered.peek(2)[:2] == b"\x1f\x8b":
                        with gzip.GzipFile(fileobj=buffered) as reader:
                            shutil.copyfileobj(reader, output, length=1024 * 1024)
                    else:
                        shutil.copyfileobj(buffered, output, length=1024 * 1024)
                output.close()
                if expected_size is not None and counted.bytes_read != expected_size:
                    raise ValueError(
                        "Genome FASTA download size differs from upstream identity metadata"
                    )
                with open(temporary, "rb") as check:
                    if check.read(1) != b">":
                        raise ValueError(
                            "Genome FASTA must start with a FASTA header: %s"
                            % self.source
                        )
                os.replace(temporary, self.materialized_path)
            finally:
                if os.path.exists(temporary):
                    os.unlink(temporary)

    def prepare(self, download=False, overwrite=False):
        """Resolve local data; only an explicit download may access the network."""
        if self.remote:
            if overwrite or self.installed_path is None:
                if not download:
                    raise MissingGenomeFastaError(
                        "Genome FASTA is not installed: %s. Call download_genome_fasta() "
                        "or run pyensembl install --only-genome-fasta." % self.source
                    )
                self.close()
                with urlopen(self.source, timeout=3600) as stream:
                    self._materialize(stream)
        else:
            if not os.path.isfile(self.source):
                raise MissingGenomeFastaError(
                    "Missing local genome FASTA: %s" % self.source
                )
            if self._local_is_compressed() and (
                overwrite or self.installed_path is None
            ):
                self.close()
                fingerprint = _fingerprint(self.source)
                with open(self.source, "rb") as stream:
                    self._materialize(stream)
                if fingerprint != _fingerprint(self.source):
                    raise ValueError(
                        "Genome FASTA changed during decompression: %s" % self.source
                    )
                _write_json(self.source_state_path, fingerprint)
        path = self.installed_path
        if path is None:
            raise MissingGenomeFastaError(
                "Empty or missing genome FASTA: %s" % self.source
            )
        return path

    def remember(self):
        """Record the attachment for CLI inspection, without changing defaults."""
        _write_json(self.manifest_path, {"source": self.source})

    def open(self, overwrite=False):
        from pyfaidx import Fasta

        path = self.prepare()
        fingerprint = _fingerprint(path)
        if self._reader is not None:
            if not overwrite and self._reader_fingerprint == fingerprint:
                return self._reader
            self.close()
        self.directory.mkdir(parents=True, exist_ok=True)
        if (
            overwrite
            or not self.index_path.exists()
            or _read_json(self.index_state_path) != fingerprint
        ):
            # Publish a complete index, never a half-written one. Keep the old
            # index available if parsing the new input fails.
            with tempfile.NamedTemporaryFile(
                dir=self.directory, delete=False
            ) as handle:
                temporary = handle.name
            os.unlink(temporary)
            try:
                with Fasta(path, indexname=temporary, strict_bounds=True) as fasta:
                    if not len(fasta.keys()):
                        raise ValueError(
                            "Genome FASTA contains no sequences: %s" % self.source
                        )
                if fingerprint != _fingerprint(path):
                    raise ValueError(
                        "Genome FASTA changed during indexing: %s" % self.source
                    )
                os.replace(temporary, self.index_path)
                _write_json(self.index_state_path, fingerprint)
            finally:
                if os.path.exists(temporary):
                    os.unlink(temporary)
        self._reader = Fasta(
            path,
            indexname=str(self.index_path),
            strict_bounds=True,
            build_index=False,
            rebuild=False,
        )
        self._reader_fingerprint = fingerprint
        return self._reader

    def close(self):
        if self._reader is not None:
            self._reader.close()
            self._reader = None
            self._reader_fingerprint = None

    @classmethod
    def installed_source(cls, cache_directory):
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
        indexed = self.index_path.exists() and _read_json(
            self.index_state_path
        ) == _fingerprint(path)
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
