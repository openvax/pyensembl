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

from os import remove
from os.path import dirname, exists, abspath, split, join

import datacache
import logging
from collections import Counter
import pickle
from .common import load_pickle, dump_pickle
from .fasta import _split_ens_version, parse_fasta_dictionary


logger = logging.getLogger(__name__)


# Bump the FASTA-pickle filename suffix when the on-disk layout of
# `_fasta_dictionary` changes so stale caches get rebuilt instead of
# silently loaded under the new code path. v2 keys versioned ENS IDs.
FASTA_PICKLE_SCHEMA_VERSION = "v2"


def lookup_sequence_with_version_fallback(sequence_data, identifier):
    """
    Look up ``identifier`` in ``sequence_data``, tolerating ENS ``.N``
    version-suffix mismatches in either direction.

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
    """
    if not identifier:
        return None
    sequence = sequence_data.get(identifier)
    if sequence is not None:
        return sequence
    if identifier.startswith("ENS") and "." in identifier:
        bare, version = _split_ens_version(identifier)
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
    Container for reference nucleotide and amino acid sequenes.
    """

    def __init__(self, fasta_paths, cache_directory_path=None):
        if type(fasta_paths) is str:
            fasta_paths = [fasta_paths]

        self.fasta_paths = [abspath(path) for path in fasta_paths]
        self.fasta_directory_paths = [split(path)[0] for path in self.fasta_paths]
        self.fasta_filenames = [split(path)[1] for path in self.fasta_paths]
        if cache_directory_path:
            self.cache_directory_paths = [cache_directory_path] * len(self.fasta_paths)
        else:
            self.cache_directory_paths = self.fasta_directory_paths
        for path in self.fasta_paths:
            if not exists(path):
                raise ValueError("Couldn't find FASTA file %s" % (path,))
        self.fasta_dictionary_filenames = [
            "%s.%s.pickle" % (filename, FASTA_PICKLE_SCHEMA_VERSION)
            for filename in self.fasta_filenames
        ]
        self.fasta_dictionary_pickle_paths = [
            join(cache_path, filename)
            for cache_path, filename in zip(
                self.cache_directory_paths, self.fasta_dictionary_filenames
            )
        ]
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

    def clear_cache(self):
        self._init_lazy_fields()
        for path in self.fasta_dictionary_pickle_paths:
            if exists(path):
                remove(path)

    def __str__(self):
        return "SequenceData(fasta_paths=%s)" % (self.fasta_paths,)

    def __repr__(self):
        return str(self)

    def __contains__(self, sequence_id):
        if self._fasta_keys is None:
            self._fasta_keys = set(self.fasta_dictionary.keys())
        return sequence_id in self._fasta_keys

    def __eq__(self, other):
        # test to see if self.fasta_paths and other.fasta_paths contain
        # the same list of paths, regardless of order
        return (other.__class__ is SequenceData) and Counter(
            self.fasta_paths
        ) == Counter(other.fasta_paths)

    def __hash__(self):
        return hash(self.fasta_paths)

    def _add_to_fasta_dictionary(self, fasta_dictionary_tmp):
        for identifier, sequence in fasta_dictionary_tmp.items():
            if identifier in self._fasta_dictionary:
                logger.warn(
                    "Sequence identifier %s is duplicated in your FASTA files!"
                    % identifier
                )
                continue
            self._fasta_dictionary[identifier] = sequence
            bare, version = _split_ens_version(identifier)
            if version is not None:
                if bare in self._stripped_index:
                    # Two FASTA entries collide on the same bare ENS ID with
                    # different versions. The newer of the two is kept as
                    # the canonical alias since version numbers grow
                    # monotonically.
                    existing = self._stripped_index[bare]
                    _, existing_version = _split_ens_version(existing)
                    if existing_version is None or version > existing_version:
                        self._stripped_index[bare] = identifier
                else:
                    self._stripped_index[bare] = identifier
                self._versions[identifier] = version

    def _load_or_create_fasta_dictionary_pickle(self):
        self._fasta_dictionary = dict()
        self._stripped_index = dict()
        self._versions = dict()
        for fasta_path, pickle_path in zip(
            self.fasta_paths, self.fasta_dictionary_pickle_paths
        ):
            if exists(pickle_path):
                # try loading the cached file
                # but we'll fall back on recreating it if loading fails
                try:
                    fasta_dictionary_tmp = load_pickle(pickle_path)
                    self._add_to_fasta_dictionary(fasta_dictionary_tmp)
                    logger.info("Loaded sequence dictionary from %s", pickle_path)
                    continue
                except (pickle.UnpicklingError, AttributeError):
                    # catch either an UnpicklingError or an AttributeError
                    # resulting from pickled objects refering to classes
                    # that no longer exists
                    logger.warn(
                        "Failed to load %s, attempting to read FASTA directly",
                        pickle_path,
                    )
            logger.info("Parsing sequences from FASTA file at %s", fasta_path)

            fasta_dictionary_tmp = parse_fasta_dictionary(fasta_path)
            self._add_to_fasta_dictionary(fasta_dictionary_tmp)
            logger.info("Saving sequence dictionary to %s", pickle_path)
            datacache.ensure_dir(dirname(pickle_path))
            dump_pickle(fasta_dictionary_tmp, pickle_path)

    def index(self, overwrite=False):
        if overwrite:
            self.clear_cache()
        self._load_or_create_fasta_dictionary_pickle()

    @property
    def fasta_dictionary(self):
        if not self._fasta_dictionary:
            self._load_or_create_fasta_dictionary_pickle()
        return self._fasta_dictionary

    def get(self, sequence_id):
        """Get sequence associated with given ID or return None if missing"""
        return self.fasta_dictionary.get(sequence_id)

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
