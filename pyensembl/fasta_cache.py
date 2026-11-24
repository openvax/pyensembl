# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Derived FASTA sequence-dictionary cache backed by ``datacache``.

Newly created FASTA sequence dictionaries are stored in a safe, versioned
gzipped JSON format (``.json.gz``) instead of Python pickles. Existing legacy
``.pickle`` caches remain readable offline and in read-only directories, and
are cleanly superseded by ``.json.gz`` files whenever an index is rebuilt.
"""

from functools import lru_cache
import gzip
import hashlib
import json
from os import lstat, remove
from os.path import abspath, dirname, join, lexists
from pathlib import Path
import pickle
import zlib

import datacache

from .common import atomic_output, file_digest, file_fingerprint
from .fasta import first_fasta_id
from .versioned_ids import split_versioned_id

FASTA_INDEX_FORMAT = "pyensembl-fasta-dictionary"
FASTA_INDEX_VERSION = 1
FASTA_INDEX_SUFFIX = ".json.gz"
LEGACY_FASTA_INDEX_SUFFIX = ".pickle"


def is_fasta_index_filename(filename):
    """Return whether ``filename`` has a current or legacy FASTA index suffix.

    Parameters
    ----------
    filename : str
        Filename or path string to check.

    Returns
    -------
    bool
        ``True`` if ``filename`` ends with :data:`FASTA_INDEX_SUFFIX` or
        :data:`LEGACY_FASTA_INDEX_SUFFIX`.
    """
    return filename.endswith((FASTA_INDEX_SUFFIX, LEGACY_FASTA_INDEX_SUFFIX))


def local_index_prefix(source):
    """Return the source-keyed filename prefix for a local FASTA cache entry.

    Parameters
    ----------
    source : str
        Local filesystem path of the source FASTA.

    Returns
    -------
    str
        Prefix of the form ``"local-<32-hex>-"``.
    """
    digest = hashlib.sha256(abspath(source).encode("utf-8")).hexdigest()[:32]
    return "local-" + digest + "-"


def legacy_local_index_prefix(source):
    """Return the pre-2.41 full-SHA-256 filename prefix for a local pickle cache.

    Parameters
    ----------
    source : str
        Local filesystem path of the source FASTA.

    Returns
    -------
    str
        Prefix of the form ``"local-<64-hex>-"``.
    """
    digest = hashlib.sha256(abspath(source).encode("utf-8")).hexdigest()
    return "local-" + digest + "-"


def safe_filename_tail(source, max_bytes=48):
    """Return a UTF-8-safe trailing slice of ``source``'s basename.

    Parameters
    ----------
    source : str
        Source path or filename.
    max_bytes : int, optional
        Maximum UTF-8 byte length of the returned suffix (default 48).

    Returns
    -------
    str
        Trailing portion of the basename fitting within ``max_bytes``.
    """
    return Path(source).name.encode("utf-8")[-max_bytes:].decode("utf-8", errors="ignore")


@lru_cache(maxsize=128)
def cached_file_digest(path, state):
    """Compute and memoize a file's SHA-256 digest while its stat identity agrees.

    Parameters
    ----------
    path : str
        Absolute path to the file.
    state : tuple
        Tuple returned by ``tuple(file_fingerprint(path, include_dev=True))``.

    Returns
    -------
    str
        Lowercase 64-character SHA-256 hex digest.

    Raises
    ------
    ValueError
        If the file's stat identity changes while hashing.
    """
    digest = file_digest(path)
    if tuple(file_fingerprint(path, include_dev=True)) != state:
        raise ValueError("FASTA changed while hashing: %s" % path)
    return digest


def file_content_digest(path):
    """Return the SHA-256 hex digest of ``path``, reusing cached values safely.

    Parameters
    ----------
    path : str
        Path to the FASTA file.

    Returns
    -------
    str
        Lowercase 64-character SHA-256 hex digest.
    """
    resolved = abspath(path)
    state = tuple(file_fingerprint(resolved, include_dev=True))
    return cached_file_digest(resolved, state)


def validate_sequence_dictionary(data, path):
    """Validate that ``data`` is a ``dict[str, str]`` sequence mapping.

    Parameters
    ----------
    data : object
        Deserialized cache payload.
    path : str or os.PathLike
        Path of the cache file (for error reporting).

    Returns
    -------
    dict of str to str
        Validated sequence dictionary.

    Raises
    ------
    datacache.FileValidationError
        If ``data`` is not a mapping of string IDs to string sequences.
    """
    if not isinstance(data, dict) or not all(
        isinstance(key, str) and isinstance(value, str)
        for key, value in data.items()
    ):
        raise datacache.FileValidationError(
            str(path), "invalid FASTA dictionary cache contents"
        )
    return data


def load_fasta_dictionary(path):
    """Load a cached FASTA sequence dictionary from ``.json.gz`` or legacy ``.pickle``.

    Parameters
    ----------
    path : str or os.PathLike
        Path to a cached FASTA dictionary file.

    Returns
    -------
    dict of str to str
        Mapping from FASTA record identifier to sequence string.

    Raises
    ------
    datacache.FileValidationError
        If the file cannot be decoded or fails format/type validation.
    """
    path_str = str(path)
    if path_str.endswith(LEGACY_FASTA_INDEX_SUFFIX):
        try:
            with open(path_str, "rb") as handle:
                data = pickle.load(handle)
        except (OSError, EOFError, ValueError, TypeError, AttributeError, pickle.UnpicklingError) as error:
            raise datacache.FileValidationError(
                path_str, "failed to read legacy pickle cache: %s" % error
            ) from error
        return validate_sequence_dictionary(data, path_str)

    try:
        with gzip.open(path_str, "rt", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, EOFError, ValueError, UnicodeDecodeError, zlib.error) as error:
        raise datacache.FileValidationError(
            path_str, "failed to read FASTA dictionary cache: %s" % error
        ) from error

    if (
        not isinstance(payload, dict)
        or payload.get("format") != FASTA_INDEX_FORMAT
        or payload.get("version") != FASTA_INDEX_VERSION
        or "sequences" not in payload
    ):
        raise datacache.FileValidationError(
            path_str, "unrecognized FASTA dictionary cache format"
        )
    return validate_sequence_dictionary(payload["sequences"], path_str)


def save_fasta_dictionary(sequences, path):
    """Atomically write ``sequences`` as a versioned gzipped JSON cache file.

    Parameters
    ----------
    sequences : dict of str to str
        Mapping from FASTA record identifier to sequence string.
    path : str or os.PathLike
        Destination path (typically ending in :data:`FASTA_INDEX_SUFFIX`).

    Returns
    -------
    str
        Destination path string written on disk.
    """
    path_str = str(path)
    validate_sequence_dictionary(sequences, path_str)
    datacache.ensure_dir(dirname(abspath(path_str)))
    payload = {
        "format": FASTA_INDEX_FORMAT,
        "version": FASTA_INDEX_VERSION,
        "sequences": sequences,
    }
    encoded = json.dumps(
        payload, ensure_ascii=True, separators=(",", ":")
    ).encode("utf-8")
    with atomic_output(path_str, "wb") as handle:
        with gzip.GzipFile(fileobj=handle, mode="wb", compresslevel=1, mtime=0) as gz:
            gz.write(encoded)
    return path_str


def cached_dictionary_lost_versions(fasta_path, fasta_dictionary):
    """Detect whether a cached FASTA dictionary dropped Ensembl header versions.

    PyEnsembl releases prior to 2.10 stripped ``.N`` version suffixes from
    Ensembl FASTA headers when building pickles. If the first FASTA header
    carries a versioned ID that is absent from ``fasta_dictionary`` while its
    bare ID is present, the cache is stale and must be rebuilt.

    Parameters
    ----------
    fasta_path : str
        Path to the source FASTA file.
    fasta_dictionary : dict of str to str
        Deserialized sequence dictionary.

    Returns
    -------
    bool
        ``True`` if the cache is missing version suffixes present in the FASTA.
    """
    try:
        first_id = first_fasta_id(fasta_path)
    except (OSError, EOFError, ValueError, zlib.error):
        # An unreadable header cannot prove the cache is stale; keep using it.
        return False
    bare, version = split_versioned_id(first_id)
    return (
        version is not None
        and first_id not in fasta_dictionary
        and bare in fasta_dictionary
    )


class FastaDictionaryCache:
    """Derived FASTA sequence-dictionary cache backed by ``datacache.Cache``.

    Centralizes cache-path naming, read-only inspection, safe gzipped-JSON
    serialization, legacy ``.pickle`` compatibility, and index cleanup under a
    single ``datacache.Cache`` root.

    Parameters
    ----------
    cache_directory : str or os.PathLike
        Directory holding derived FASTA dictionary indexes.
    """

    def __init__(self, cache_directory):
        self.cache = datacache.Cache(cache_root=str(cache_directory))

    @property
    def cache_directory(self):
        """Resolved cache root directory string."""
        return self.cache.cache_directory_path

    def target_path(self, fasta_path, source=None, *, content_path=None):
        """Return the canonical ``.json.gz`` cache path for a FASTA source.

        Never touches the filesystem except to hash a local FASTA's bytes so
        local edits select a distinct cache generation.

        Parameters
        ----------
        fasta_path : str
            Local filesystem path to the FASTA file.
        source : str, optional
            Original configured path or URL (defaults to ``fasta_path``).
        content_path : str, optional
            Path whose bytes determine the content digest for local sources
            (defaults to ``fasta_path``).

        Returns
        -------
        str
            Destination ``.json.gz`` path inside :attr:`cache_directory`.
        """
        source = source or fasta_path
        if "://" in source:
            base = Path(fasta_path).name
            candidate = base + FASTA_INDEX_SUFFIX
            if len(candidate) > 150 or len(candidate.encode("utf-8")) > 180:
                digest = hashlib.sha256(base.encode("utf-8")).hexdigest()[:32]
                candidate = (
                    "remote-" + digest + "-" + safe_filename_tail(base, 48) + FASTA_INDEX_SUFFIX
                )
            return self.cache.local_path(filename=candidate)

        content = abspath(content_path or fasta_path)
        try:
            digest = file_content_digest(content)[:32]
        except OSError:
            # Read-only inspection of a missing/unreadable source stays side-effect-free.
            digest = "missing"
        candidate = (
            local_index_prefix(source)
            + digest
            + "-"
            + safe_filename_tail(source, 48)
            + FASTA_INDEX_SUFFIX
        )
        return self.cache.local_path(filename=candidate)

    def legacy_path(self, fasta_path, source=None, *, content_path=None):
        """Locate an existing legacy ``.pickle`` cache for ``fasta_path``, if any.

        Parameters
        ----------
        fasta_path : str
            Local filesystem path to the FASTA file.
        source : str, optional
            Original configured path or URL (defaults to ``fasta_path``).
        content_path : str, optional
            Path whose bytes determine the content digest for local sources.

        Returns
        -------
        str or None
            Existing legacy ``.pickle`` path, or ``None`` if none applies.
        """
        source = source or fasta_path
        directory = Path(self.cache_directory)
        legacy_basename = join(
            self.cache_directory, Path(fasta_path).name + LEGACY_FASTA_INDEX_SUFFIX
        )
        if "://" in source:
            return legacy_basename if lexists(legacy_basename) else None

        content = abspath(content_path or fasta_path)
        try:
            full_digest = file_content_digest(content)
        except OSError:
            full_digest = None
        if full_digest is not None:
            old_tail = safe_filename_tail(source, 72)
            hashed_legacy = join(
                self.cache_directory,
                legacy_local_index_prefix(source)
                + full_digest
                + "-"
                + old_tail
                + LEGACY_FASTA_INDEX_SUFFIX,
            )
            if lexists(hashed_legacy):
                return hashed_legacy

        has_keyed = any(
            directory.glob(local_index_prefix(source) + "*" + FASTA_INDEX_SUFFIX)
        ) or any(
            directory.glob(legacy_local_index_prefix(source) + "*" + LEGACY_FASTA_INDEX_SUFFIX)
        )
        if lexists(legacy_basename) and not has_keyed:
            return legacy_basename
        return None

    def index_path(self, fasta_path, source=None, *, content_path=None):
        """Select the active index path without creating or migrating files.

        Prefers an existing ``.json.gz`` index, falls back to an existing
        legacy ``.pickle`` index when present, and otherwise returns the
        canonical ``.json.gz`` target path.

        Parameters
        ----------
        fasta_path : str
            Local filesystem path to the FASTA file.
        source : str, optional
            Original configured path or URL (defaults to ``fasta_path``).
        content_path : str, optional
            Path whose bytes determine the content digest for local sources.

        Returns
        -------
        str
            Selected cache file path inside :attr:`cache_directory`.
        """
        target = self.target_path(fasta_path, source, content_path=content_path)
        if lexists(target):
            return target
        legacy = self.legacy_path(fasta_path, source, content_path=content_path)
        if legacy is not None:
            return legacy
        return target

    def index_files(self, fasta_path, source=None):
        """Return all existing index files owned by ``source`` (current and legacy).

        Includes older byte-digest generations for local sources so index
        deletion removes every derived cache owned by this source.

        Parameters
        ----------
        fasta_path : str
            Local filesystem path to the FASTA file.
        source : str, optional
            Original configured path or URL (defaults to ``fasta_path``).

        Returns
        -------
        list of str
            Existing index file paths owned by this source.
        """
        directory = Path(self.cache_directory)
        source = source or fasta_path
        if "://" not in source:
            paths = sorted(
                list(directory.glob(local_index_prefix(source) + "*" + FASTA_INDEX_SUFFIX))
                + list(directory.glob(legacy_local_index_prefix(source) + "*" + LEGACY_FASTA_INDEX_SUFFIX))
            )
            if paths:
                return [str(path) for path in paths]
            legacy = directory / (Path(fasta_path).name + LEGACY_FASTA_INDEX_SUFFIX)
            return [str(legacy)] if lexists(legacy) else []

        candidates = [
            self.target_path(fasta_path, source),
            join(self.cache_directory, Path(fasta_path).name + LEGACY_FASTA_INDEX_SUFFIX),
        ]
        return [path for path in dict.fromkeys(candidates) if lexists(path)]

    def inspect(self, fasta_path, source=None, *, content_path=None):
        """Inspect the selected FASTA index file offline via ``datacache.inspect_file``.

        Parameters
        ----------
        fasta_path : str
            Local filesystem path to the FASTA file.
        source : str, optional
            Original configured path or URL (defaults to ``fasta_path``).
        content_path : str, optional
            Path whose bytes determine the content digest for local sources.

        Returns
        -------
        datacache.FileInspection
            Read-only file inspection result.
        """
        return datacache.inspect_file(
            self.index_path(fasta_path, source, content_path=content_path)
        )

    def load(self, path):
        """Load a cached FASTA dictionary from ``path``.

        Parameters
        ----------
        path : str or os.PathLike
            Path to a ``.json.gz`` or legacy ``.pickle`` FASTA dictionary.

        Returns
        -------
        dict of str to str
            Mapping from sequence identifier to sequence string.
        """
        return load_fasta_dictionary(path)

    def save(self, sequences, fasta_path, source=None, *, content_path=None):
        """Save ``sequences`` to the canonical ``.json.gz`` cache path.

        Parameters
        ----------
        sequences : dict of str to str
            Mapping from sequence identifier to sequence string.
        fasta_path : str
            Local filesystem path to the FASTA file.
        source : str, optional
            Original configured path or URL (defaults to ``fasta_path``).
        content_path : str, optional
            Path whose bytes determine the content digest for local sources.

        Returns
        -------
        str
            Path string of the published ``.json.gz`` cache file.
        """
        target = self.target_path(fasta_path, source, content_path=content_path)
        return save_fasta_dictionary(sequences, target)

    def delete(self, fasta_path, source=None):
        """Delete all existing index files owned by ``source``.

        Parameters
        ----------
        fasta_path : str
            Local filesystem path to the FASTA file.
        source : str, optional
            Original configured path or URL (defaults to ``fasta_path``).

        Returns
        -------
        list of tuple[str, int]
            ``(path, size_in_bytes)`` pairs for each removed index file.
        """
        deleted = []
        for path in self.index_files(fasta_path, source):
            if lexists(path):
                size = lstat(path).st_size
                remove(path)
                deleted.append((path, size))
        return deleted


def fasta_index_files(fasta_path, cache_directory=None, source=None):
    """Return existing FASTA index files owned by ``source``, including older revisions.

    Parameters
    ----------
    fasta_path : str
        Local filesystem path to the FASTA file.
    cache_directory : str, optional
        Cache directory containing the indexes (defaults to ``dirname(fasta_path)``).
    source : str, optional
        Original configured path or URL (defaults to ``fasta_path``).

    Returns
    -------
    list of str
        Existing index paths (both ``.json.gz`` and legacy ``.pickle``).
    """
    cache = FastaDictionaryCache(cache_directory or dirname(fasta_path))
    return cache.index_files(fasta_path, source)


def fasta_index_path(fasta_path, cache_directory=None, source=None, *, content_path=None):
    """Select the active FASTA index path without creating or migrating files.

    Parameters
    ----------
    fasta_path : str
        Local filesystem path to the FASTA file.
    cache_directory : str, optional
        Cache directory containing the index (defaults to ``dirname(fasta_path)``).
    source : str, optional
        Original configured path or URL (defaults to ``fasta_path``).
    content_path : str, optional
        Path whose bytes determine the content digest for local sources.

    Returns
    -------
    str
        Active index path inside the cache directory.
    """
    cache = FastaDictionaryCache(cache_directory or dirname(fasta_path))
    return cache.index_path(fasta_path, source, content_path=content_path)
