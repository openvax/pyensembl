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

from collections import Counter
import logging
from os.path import abspath, exists, split
from pathlib import Path

from datacache import FileValidationError

from .fasta import parse_fasta_dictionary
from .fasta_cache import (  # noqa: F401
    FastaDictionaryCache,
    cached_dictionary_lost_versions,
    fasta_index_files,
    fasta_index_path,
)
from .versioned_ids import match_version, split_versioned_id


logger = logging.getLogger(__name__)


def lookup_sequence_with_version_fallback(sequence_data, identifier, version=None):
    """
    Look up ``identifier`` in ``sequence_data``, tolerating ENS ``.N``
    version-suffix mismatches in either direction. ``version``, the version
    the GTF records for a bare ``identifier``, is preferred when the FASTA
    has that record.

    Both Ensembl and GENCODE work with versioned IDs (e.g.
    ``ENSP00000123456.3``); the two formats just split that information
    differently. Ensembl GTFs keep the bare ID in ``protein_id`` /
    ``transcript_id`` and put the version in a separate ``*_version``
    attribute, while GENCODE GTFs embed the version directly in the ID
    attribute. Whichever convention the FASTA on disk used, the
    sequence dict now keys on the form that header carried — so either
    a versioned or a bare caller-supplied ID may need a fallback.

    Resolution order:
      1. literal lookup of ``identifier``
      2. if ``identifier`` is a versioned ENS ID, strip ``.N`` and retry
         (handles caller-supplied versioned ID against a bare FASTA)
      3. consult ``sequence_data._stripped_index`` for a versioned alias
         of a bare caller-supplied ID (handles bare ID against a
         versioned FASTA, the GENCODE case)
      4. None

    Non-Ensembl IDs (e.g. TAIR ``AT1G01010.1`` where ``.N`` is an
    isoform suffix, not a version) skip step (2): stripping there is
    semantically wrong.

    Use this for IDs the annotation itself stores (object properties such
    as ``Transcript.sequence``). IDs from callers go through
    :meth:`SequenceData.stored_id`, which checks a supplied version.
    """
    if not identifier:
        return None
    if version is not None:
        sequence = sequence_data.get("%s.%d" % (identifier, version))
        if sequence is not None:
            return sequence
    sequence = sequence_data.get(identifier)
    if sequence is not None:
        return sequence
    if identifier.startswith("ENS") and "." in identifier:
        bare, version = split_versioned_id(identifier)
        if version is not None:
            sequence = sequence_data.get(bare)
            if sequence is not None:
                return sequence
    stripped_index = getattr(sequence_data, "_stripped_index", None)
    if stripped_index:
        versioned = stripped_index.get(identifier)
        if versioned is not None:
            return sequence_data.get(versioned)
    return None


class SequenceData(object):
    """
    Container for reference nucleotide and amino acid sequences.
    """

    def __init__(self, fasta_paths, cache_directory_path=None, *, sources=None):
        if type(fasta_paths) is str:
            fasta_paths = [fasta_paths]

        self.fasta_paths = [abspath(path) for path in fasta_paths]
        self.sources = list(sources) if sources is not None else self.fasta_paths
        if len(self.sources) != len(self.fasta_paths):
            raise ValueError("Expected one source for each FASTA path")
        self.fasta_directory_paths = [split(path)[0] for path in self.fasta_paths]
        self.fasta_filenames = [split(path)[1] for path in self.fasta_paths]
        if cache_directory_path:
            self.cache_directory_paths = [cache_directory_path] * len(self.fasta_paths)
        else:
            self.cache_directory_paths = self.fasta_directory_paths
        for path in self.fasta_paths:
            if not exists(path):
                raise ValueError("Couldn't find FASTA file %s" % (path,))
        self.fasta_dictionary_index_paths = []
        self._refresh_index_paths()
        self._init_lazy_fields()

    @property
    def fasta_dictionary_pickle_paths(self):
        """Backwards-compatible alias for :attr:`fasta_dictionary_index_paths`."""
        return self.fasta_dictionary_index_paths

    def _refresh_index_paths(self):
        paths = [
            fasta_index_path(path, directory, source)
            for path, directory, source in zip(
                self.fasta_paths,
                self.cache_directory_paths,
                getattr(self, "sources", self.fasta_paths),
            )
        ]
        if paths != self.fasta_dictionary_index_paths:
            self.fasta_dictionary_index_paths = paths
            self.fasta_dictionary_filenames = [Path(path).name for path in paths]
            self._init_lazy_fields()

    def _init_lazy_fields(self):
        self._fasta_dictionary = None
        self._fasta_keys = None
        # Maps a bare ENS ID -> the versioned form actually keyed in
        # _fasta_dictionary, so callers passing the unversioned form can
        # still resolve to a versioned FASTA entry (GENCODE case).
        self._stripped_index = None
        # Maps a sequence identifier -> the integer version suffix the
        # FASTA header carried, or absent for headers without a version.
        # Lets callers sanity-check FASTA / GTF alignment.
        self._versions = None
        # Maps a bare ENS ID with several versioned FASTA entries -> all of
        # them, so a bare lookup can report the ambiguity
        self._version_conflicts = None

    def clear_cache(self):
        """Clear parsed sequence data held in memory."""
        self._init_lazy_fields()

    def delete_index_files(self):
        """Delete cached FASTA dictionaries while preserving source files."""
        self.clear_cache()
        for fasta_path, directory, source in zip(
            self.fasta_paths,
            self.cache_directory_paths,
            getattr(self, "sources", self.fasta_paths),
        ):
            FastaDictionaryCache(directory).delete(fasta_path, source)
        self._refresh_index_paths()

    def __str__(self):
        return "SequenceData(fasta_paths=%s)" % (self.fasta_paths,)

    def __repr__(self):
        return str(self)

    def __contains__(self, sequence_id):
        _ = self.fasta_dictionary
        if self._fasta_keys is None:
            self._fasta_keys = set(self.fasta_dictionary.keys())
        return sequence_id in self._fasta_keys

    def __eq__(self, other):
        # test to see if self.fasta_paths and other.fasta_paths contain
        # the same list of paths, regardless of order
        return (type(self) is type(other)) and Counter(
            self.fasta_paths
        ) == Counter(other.fasta_paths)

    def __hash__(self):
        # Match equality: ignore path order but retain duplicate counts.
        return hash(tuple(sorted(self.fasta_paths)))

    def _add_to_fasta_dictionary(self, fasta_dictionary_tmp):
        for identifier, sequence in fasta_dictionary_tmp.items():
            if identifier in self._fasta_dictionary:
                logger.warning(
                    "Sequence identifier %s is duplicated in your FASTA files!"
                    % identifier
                )
                continue
            self._fasta_dictionary[identifier] = sequence
            bare, version = split_versioned_id(identifier)
            if version is not None:
                if bare in self._stripped_index:
                    # Two FASTA entries collide on the same bare ENS ID with
                    # different versions. The newer of the two is kept as
                    # the canonical alias since version numbers grow
                    # monotonically.
                    existing = self._stripped_index[bare]
                    self._version_conflicts.setdefault(bare, {existing}).add(
                        identifier
                    )
                    _, existing_version = split_versioned_id(existing)
                    if existing_version is None or version > existing_version:
                        self._stripped_index[bare] = identifier
                else:
                    self._stripped_index[bare] = identifier
                self._versions[identifier] = version

    def load_or_create_fasta_dictionary(self, show_progress=False):
        """Load cached FASTA sequence dictionaries or parse and cache them.

        Reads existing ``.json.gz`` or legacy ``.pickle`` caches through
        :class:`~pyensembl.fasta_cache.FastaDictionaryCache`. Newly built or
        rebuilt caches are always published in the safe ``.json.gz`` format.

        Parameters
        ----------
        show_progress : bool, optional
            Display a progress bar while parsing FASTA files from disk.
        """
        self._refresh_index_paths()
        self._fasta_dictionary = dict()
        self._stripped_index = dict()
        self._versions = dict()
        self._version_conflicts = dict()
        for idx, (fasta_path, index_path, directory, source) in enumerate(
            zip(
                self.fasta_paths,
                self.fasta_dictionary_index_paths,
                self.cache_directory_paths,
                getattr(self, "sources", self.fasta_paths),
            )
        ):
            cache = FastaDictionaryCache(directory)
            if exists(index_path):
                try:
                    fasta_dictionary_tmp = cache.load(index_path)
                    if cached_dictionary_lost_versions(fasta_path, fasta_dictionary_tmp):
                        logger.info(
                            "Rebuilding %s to keep FASTA header versions",
                            index_path,
                        )
                    else:
                        self._add_to_fasta_dictionary(fasta_dictionary_tmp)
                        logger.debug(
                            "Loaded sequence dictionary from %s", index_path
                        )
                        continue
                except FileValidationError:
                    logger.warning(
                        "Failed to load %s, attempting to read FASTA directly",
                        index_path,
                    )
            logger.info("Parsing sequences from FASTA file at %s", fasta_path)

            expected_target = cache.target_path(fasta_path, source)
            fasta_dictionary_tmp = parse_fasta_dictionary(fasta_path, show_progress)
            if cache.target_path(fasta_path, source) != expected_target:
                self._init_lazy_fields()
                raise ValueError("FASTA changed while indexing: %s" % fasta_path)
            self._add_to_fasta_dictionary(fasta_dictionary_tmp)
            logger.debug("Saving sequence dictionary to %s", expected_target)
            try:
                saved_path = cache.save(fasta_dictionary_tmp, fasta_path, source)
                self.fasta_dictionary_index_paths[idx] = saved_path
                self.fasta_dictionary_filenames[idx] = Path(saved_path).name
            except OSError as e:
                # A read-only shared cache still serves the parsed sequences
                logger.warning("Couldn't save %s: %s", expected_target, e)

    def _load_or_create_fasta_dictionary_pickle(self, show_progress=False):
        self.load_or_create_fasta_dictionary(show_progress=show_progress)

    def index(self, overwrite=False, show_progress=False):
        if overwrite:
            self.delete_index_files()
        self.load_or_create_fasta_dictionary(show_progress=show_progress)

    @property
    def fasta_dictionary(self):
        self._refresh_index_paths()
        if self._fasta_dictionary is None:
            self.load_or_create_fasta_dictionary()
        return self._fasta_dictionary

    def get(self, sequence_id):
        """Get sequence associated with given ID or return None if missing"""
        return self.fasta_dictionary.get(sequence_id)

    def stored_id(self, sequence_id):
        """
        FASTA record ID for ``sequence_id``, or ``None`` if it has no record.

        A bare Ensembl ID matches the record whatever its version, and raises
        ``ValueError`` if several versions are present. A versioned ID must
        match the header's version: ``ValueError`` if the header carries a
        different version or none.
        """
        if not sequence_id:
            return None
        if sequence_id in self.fasta_dictionary:
            return sequence_id
        bare, _ = split_versioned_id(sequence_id)
        installed = {form: self._versions[form] for form in self._versioned_forms(bare)}
        if bare in self.fasta_dictionary:
            installed[bare] = None
        return match_version(sequence_id, installed)

    def recorded_versions(self, sequence_id):
        """Sorted versions that FASTA headers record for this stable ID."""
        bare, _ = split_versioned_id(sequence_id)
        return sorted(self._versions[form] for form in self._versioned_forms(bare))

    def _versioned_forms(self, bare):
        _ = self.fasta_dictionary  # ensure lazy load
        if bare in self._version_conflicts:
            return self._version_conflicts[bare]
        if bare in self._stripped_index:
            return [self._stripped_index[bare]]
        return []

    def fasta_version(self, sequence_id):
        """
        Return the integer FASTA-header version for ``sequence_id``, or
        ``None`` if the header didn't carry a version (e.g. older Ensembl
        releases) or the ID isn't in the FASTA.

        Accepts either the versioned or bare form of an ENS ID — the
        bare form is resolved through ``_stripped_index``.
        """
        if not sequence_id:
            return None
        # ensure lazy load
        _ = self.fasta_dictionary
        version = self._versions.get(sequence_id)
        if version is not None:
            return version
        versioned = self._stripped_index.get(sequence_id)
        if versioned is not None:
            return self._versions.get(versioned)
        return None
