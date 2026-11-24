"""Disk-backed reference DNA, with indexes owned by PyEnsembl."""

import gzip
import hashlib
import logging
import os
from pathlib import Path
import re
from urllib.parse import urlsplit

from datacache import FileValidationError, fetch_file

from .common import (
    atomic_output,
    file_digest,
    file_fingerprint,
    publish_file,
    read_json,
    read_progress,
    remove_if_exists,
    staging_path,
    write_json,
)
from .download_cache import DOWNLOAD_TIMEOUT
from .normalization import normalize_chromosome

logger = logging.getLogger(__name__)
_CHUNK_SIZE = 1024 * 1024

_UCSC_UNPLACED_RE = re.compile(r"^chrun_([a-z0-9]+)v(\d+)$")
_UCSC_SCAFFOLD_RE = re.compile(r"^chr[a-z0-9]+_([a-z0-9]+)v(\d+)_(?:random|alt|fix)$")
_MITO_ALIASES = frozenset({"m", "mt", "mito"})


def contig_alias_key(contig):
    """Return a canonical syntactic alias key for comparing FASTA contig names.

    Recognizes common Ensembl/UCSC/GENCODE naming conventions without altering
    accession versions or assuming assembly identity:

    - Optional case-insensitive ``chr`` prefix (e.g. ``1`` <-> ``chr1``,
      ``X`` <-> ``chrX``) when followed by a non-underscore character.
    - Mitochondrial synonyms (``M``, ``MT``, ``Mito``, ``chrM``, ``chrMT``,
      ``chrMito`` -> ``"mt"``).
    - Version-preserving UCSC unplaced/unlocalized/alt/fix scaffold names
      (e.g. ``chrUn_KI270728v1`` <-> ``KI270728.1``,
      ``chr4_GL000194v1_random`` <-> ``GL000194.1``,
      ``chr1_KI270758v1_alt`` <-> ``KI270758.1``).

    Parameters
    ----------
    contig : str or int
        Contig or chromosome name.

    Returns
    -------
    str
        Lowercase alias key suitable for grouping syntactically equivalent
        contig spellings.
    """
    name = str(contig).strip().lower()
    unplaced = _UCSC_UNPLACED_RE.fullmatch(name)
    if unplaced is not None:
        return "%s.%s" % (unplaced.group(1), unplaced.group(2))
    scaffold = _UCSC_SCAFFOLD_RE.fullmatch(name)
    if scaffold is not None:
        return "%s.%s" % (scaffold.group(1), scaffold.group(2))
    if name.startswith("chr") and len(name) > 3 and name[3] != "_":
        name = name[3:]
    if name in _MITO_ALIASES:
        return "mt"
    return name


def format_contig_coverage_warning(coverage):
    """Format a user warning for genuinely missing or ambiguous FASTA contigs.

    Contigs that match unambiguously via naming aliases (e.g. ``1`` -> ``chr1``
    or ``MT`` -> ``chrM``) do not trigger a warning by themselves, while
    genuinely unavailable or ambiguously aliased contigs produce a clear
    warning that distinguishes them from alias-compatible naming differences.

    Parameters
    ----------
    coverage : dict
        Result dictionary returned by :meth:`GenomeFasta.classify_contigs` with
        keys ``"matched"``, ``"aliased"``, ``"ambiguous"``, and ``"missing"``.

    Returns
    -------
    str or None
        Warning message to emit, or ``None`` if no contigs are missing or
        ambiguous.
    """
    missing = sorted(coverage["missing"])
    ambiguous = coverage["ambiguous"]
    aliased = coverage["aliased"]
    if not missing and not ambiguous:
        return None
    parts = []
    if missing:
        parts.append(
            "Local genome FASTA lacks %d annotation contigs (e.g. %s)"
            % (len(missing), ", ".join(missing[:5]))
        )
    if ambiguous:
        ambiguous_examples = ", ".join(
            "%s (%s)" % (name, "/".join(ambiguous[name]))
            for name in sorted(ambiguous)[:5]
        )
        if missing:
            parts.append(
                "%d annotation contigs match multiple FASTA aliases ambiguously (e.g. %s)"
                % (len(ambiguous), ambiguous_examples)
            )
        else:
            parts.append(
                "Local genome FASTA has %d ambiguously aliased annotation contigs (e.g. %s)"
                % (len(ambiguous), ambiguous_examples)
            )
    if aliased:
        alias_examples = ", ".join(
            "%s->%s" % (name, aliased[name]) for name in sorted(aliased)[:3]
        )
        parts.append(
            "%d other annotation contigs match only via naming aliases (e.g. %s)"
            % (len(aliased), alias_examples)
        )
    return (
        "; ".join(parts)
        + ". Check assembly and contig naming; matching names or aliases alone "
        "do not verify assembly identity."
    )


class MissingGenomeFastaError(ValueError):
    """Reference DNA was not configured or has not been downloaded."""


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
        # A different attachment may now occupy the advisory manifest. Keep
        # reopening an existing source read-only, including pre-receipt caches.
        self._attachment_known = self.materialized_path.is_file() or self.index_path.is_file()
        if not self._attachment_known:
            state = {"source": self.source}
            manifest = read_json(self.manifest_path)
            self._attachment_known = (
                read_json(self.directory / "attachment.json") == state
                or isinstance(manifest, dict) and manifest.get("source") == self.source
            )
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
        return file_fingerprint(path, ctime=path == self.source)

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
            if read_json(self.source_state_path) != file_fingerprint(self.source):
                return None
        return path

    def _not_installed(self):
        message = "Reference DNA is not installed: %s. Call download_genome_fasta()" % (
            self.source
        )
        if self.install_string_function is not None:
            message += " or run: %s" % self.install_string_function()
        return MissingGenomeFastaError(message)

    def _write_uncompressed(self, source_path, destination, show_progress=False):
        """Publish plain FASTA from a plain or gzip file, or nothing on failure."""
        description = "Decompressing " + os.path.basename(urlsplit(self.source).path)
        with open(source_path, "rb") as raw, read_progress(
            raw, description, show_progress
        ) as report_progress:
            compressed = raw.peek(2)[:2] == b"\x1f\x8b"
            reader = gzip.GzipFile(fileobj=raw) if compressed else raw
            with reader:
                first = reader.read(_CHUNK_SIZE)
                if not first.startswith(b">"):
                    raise ValueError(
                        "Genome FASTA must start with a FASTA header: %s" % self.source
                    )
                with atomic_output(destination) as output:
                    chunk = first
                    while chunk:
                        output.write(chunk)
                        report_progress()
                        chunk = reader.read(_CHUNK_SIZE)

    def _download(self, expected_size=None, show_progress=False, *, expected_md5=None):
        """Fetch with datacache retries, then publish the decompressed FASTA.

        expected_size describes the downloaded (compressed) bytes. When it is
        known, an interrupted download resumes where it stopped, even in a
        later process: datacache appends only bytes the server confirms are
        from the same file (strong ETag), then checks the size. A dated
        platform MD5, when supplied, is checked before publishing any FASTA.
        """
        self.directory.mkdir(parents=True, exist_ok=True)
        # A stable name lets a later attempt find the partial download, and a
        # complete one that was interrupted before decompression.
        name = os.path.basename(urlsplit(self.source).path) or "download"
        raw = self.directory / ("download-" + name)
        options = dict(
            destination=raw,
            raw=True,
            timeout=DOWNLOAD_TIMEOUT,
            expected_size=expected_size,
            show_progress=show_progress,
        )
        resume = expected_size is not None and os.name == "posix"
        try:
            logger.info("Downloading genome FASTA %s", self.source)
            try:
                fetch_file(self.source, resume=resume, **options)
            except FileValidationError as error:
                if not resume or "strong ETag" not in error.reason:
                    raise
                logger.warning(
                    "Cannot resume %s (no strong ETag); downloading it in full",
                    self.source,
                )
                fetch_file(self.source, **options)
            if expected_md5 is not None:
                if file_digest(raw, "md5") != expected_md5:
                    raise FileValidationError(str(raw), "MD5 does not match upstream manifest")
            self._write_uncompressed(raw, self.materialized_path, show_progress)
        finally:
            remove_if_exists(raw)

    def prepare(self, download=False, overwrite=False, show_progress=False):
        """Resolve local data; only an explicit download may access the network."""
        if self.remote:
            if overwrite or self.installed_path is None:
                if not download:
                    raise self._not_installed()
                self._download(show_progress=show_progress)
        else:
            if not os.path.isfile(self.source):
                raise MissingGenomeFastaError(
                    "Missing local genome FASTA: %s" % self.source
                )
            if self._local_is_compressed() and (
                overwrite or self.installed_path is None
            ):
                logger.info("Decompressing genome FASTA %s", self.source)
                fingerprint = file_fingerprint(self.source)
                self._write_uncompressed(self.source, self.materialized_path, show_progress)
                if fingerprint != file_fingerprint(self.source):
                    raise ValueError(
                        "Genome FASTA changed during decompression: %s" % self.source
                    )
                write_json(self.source_state_path, fingerprint)
        path = self.installed_path
        if path is None:
            raise self._not_installed() if self.remote else MissingGenomeFastaError(
                "Empty or missing genome FASTA: %s" % self.source
            )
        return path

    def remember(self):
        """Record the attachment for CLI inspection, without changing defaults."""
        if self._attachment_known:
            return
        state = {"source": self.source}
        if read_json(self.manifest_path) != state:
            write_json(self.manifest_path, state)
        if not self.materialized_path.is_file() and not self.index_path.is_file():
            # Download-only local input has no owned FASTA/index yet. Record it
            # per source so later attachments cannot make it appear new again.
            write_json(self.directory / "attachment.json", state)
        self._attachment_known = True

    def _index_is_current(self, fingerprint):
        return self.index_path.exists() and read_json(self.index_state_path) == fingerprint

    def _ensure_index(self, path, fingerprint, overwrite=False):
        """Publish a complete index, keeping the old one if parsing fails."""
        from pyfaidx import Fasta

        if not overwrite and self._index_is_current(fingerprint):
            return
        logger.info("Indexing genome FASTA %s", path)
        self.directory.mkdir(parents=True, exist_ok=True)
        staged = staging_path(self.directory, self.index_path.name)
        try:
            with Fasta(path, indexname=str(staged), strict_bounds=True) as fasta:
                if not len(fasta.keys()):
                    raise ValueError(
                        "Genome FASTA contains no sequences: %s" % self.source
                    )
            if fingerprint != self._file_fingerprint(path):
                raise ValueError("Genome FASTA changed during indexing: %s" % self.source)
            publish_file(staged, self.index_path)
            write_json(self.index_state_path, fingerprint)
        finally:
            remove_if_exists(staged)

    def open(self, overwrite=False, show_progress=False):
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
        path = self.prepare(show_progress=show_progress)
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
        Raises ValueError if a normalized name matches several records.
        """
        fasta = self.open()
        name = str(contig)
        if name in fasta:
            return fasta[name]
        names_by_normalized = self._names_by_normalized
        if names_by_normalized is None:
            # Build fully before publishing, for concurrent readers.
            names_by_normalized = {}
            for record in fasta.keys():
                try:
                    key = normalize_chromosome(record)
                except (TypeError, ValueError):
                    continue
                names_by_normalized.setdefault(key, []).append(record)
            self._names_by_normalized = names_by_normalized
        try:
            matches = names_by_normalized.get(normalize_chromosome(contig), [])
        except (TypeError, ValueError):
            return None
        if len(matches) > 1:
            raise ValueError(
                "Contig %r matches several FASTA records: %s" % (contig, ", ".join(matches))
            )
        return fasta[matches[0]] if matches else None

    def similar_names(self, contig):
        """Return FASTA record names matching ``contig`` via syntactic aliases.

        Matches up to case, optional ``chr`` prefix, mitochondrial synonyms
        (``MT`` / ``chrM`` / ``M``), and version-preserving UCSC scaffold
        names (e.g. ``chrUn_KI270728v1`` <-> ``KI270728.1``).

        Parameters
        ----------
        contig : str or int
            Query contig name.

        Returns
        -------
        list of str
            FASTA record names sharing :func:`contig_alias_key` with ``contig``.
        """
        names_by_alias = self._get_names_by_alias()
        return list(names_by_alias.get(contig_alias_key(contig), ()))

    def _get_names_by_alias(self):
        names_by_alias = self._names_by_alias
        if names_by_alias is None:
            fasta = self.open()
            names_by_alias = {}
            for record in fasta.keys():
                names_by_alias.setdefault(contig_alias_key(record), []).append(record)
            self._names_by_alias = names_by_alias
        return names_by_alias

    def classify_contigs(self, contigs):
        """Classify annotation contigs against the records in this genome FASTA.

        Distinguishes exact/normalized matches from unambiguous naming aliases
        (such as ``1`` <-> ``chr1`` or ``MT`` <-> ``chrM``), ambiguous alias
        collisions, and genuinely missing contigs.

        Parameters
        ----------
        contigs : iterable of str
            Annotation contig names (e.g. from ``Genome.contigs()``).

        Returns
        -------
        dict
            Dictionary with keys:

            - ``"matched"`` : dict[str, str]
              Contigs resolved directly via :meth:`record` (exact or
              ``normalize_chromosome`` match).
            - ``"aliased"`` : dict[str, str]
              Contigs not matched directly that uniquely match a FASTA record
              via :func:`contig_alias_key`.
            - ``"ambiguous"`` : dict[str, list[str]]
              Contigs matching multiple FASTA records by normalized name or
              alias key.
            - ``"missing"`` : list[str]
              Sorted contigs with no matching or aliased FASTA record.
        """
        fasta = self.open()
        names_by_alias = self._get_names_by_alias()
        matched = {}
        aliased = {}
        ambiguous = {}
        missing = []
        for raw_contig in contigs:
            contig = str(raw_contig)
            if contig in fasta:
                matched[contig] = contig
                continue
            try:
                record = self.record(contig)
            except ValueError:
                normalized = normalize_chromosome(contig)
                candidates = list(self._names_by_normalized.get(normalized, ()))
                ambiguous[contig] = candidates
                continue
            if record is not None:
                matched[contig] = record.name
                continue
            candidates = list(names_by_alias.get(contig_alias_key(contig), ()))
            if len(candidates) == 1:
                aliased[contig] = candidates[0]
            elif len(candidates) > 1:
                ambiguous[contig] = candidates
            else:
                missing.append(contig)
        return {
            "matched": matched,
            "aliased": aliased,
            "ambiguous": ambiguous,
            "missing": sorted(missing),
        }

    def _forget_reader(self):
        self._reader = None
        self._reader_paths = ()
        self._reader_fingerprints = None
        self._names_by_normalized = None
        self._names_by_alias = None

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
        state = read_json(Path(cache_directory) / "genome_fasta.json")
        if isinstance(state, dict) and isinstance(state.get("source"), str):
            if "shared_key" in state:
                from .genome_fasta_cache import SharedGenomeFasta

                return SharedGenomeFasta(state["source"], cache_directory)
            return cls(state["source"], cache_directory)
        return None

    def status(self, check=False):
        """'indexed', 'needs index', 'invalid index' (check=True), or 'missing'."""
        path = self.installed_path
        if path is None:
            return "missing"
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
                return "invalid index"
        return "indexed" if indexed else "needs index"
