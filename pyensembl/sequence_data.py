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
from .fasta import parse_fasta_dictionary


logger = logging.getLogger(__name__)


def lookup_sequence_with_version_fallback(sequence_data, identifier):
    """
    Look up ``identifier`` in ``sequence_data``, tolerating an ENS.N version
    suffix mismatch between the caller's ID and the FASTA's dict keys.

    Both Ensembl and GENCODE work with versioned IDs (e.g.
    ``ENSP00000123456.3``); the two formats just split that information
    differently. Ensembl GTFs keep the bare ID in ``protein_id`` /
    ``transcript_id`` and put the version in a separate ``*_version``
    attribute, while GENCODE GTFs embed the version directly in the ID
    attribute. pyensembl's FASTA parser strips ``.N`` from ENS-prefix
    headers, so the FASTA dict keys are always unversioned — meaning a
    literal lookup of a GENCODE-style versioned ID would miss.

    Strategy: try ``identifier`` as-is first (covers both unversioned IDs
    and any future FASTA that preserves versions); on miss, strip a
    trailing ``.N`` suffix and try again, but only for ENS-prefix IDs
    because for non-Ensembl IDs (e.g. TAIR ``AT1G01010.1``) the ``.N`` is
    an isoform suffix and stripping would be incorrect.
    """
    if not identifier:
        return None
    sequence = sequence_data.get(identifier)
    if sequence is not None:
        return sequence
    if identifier.startswith("ENS") and "." in identifier:
        stripped = identifier.rsplit(".", 1)[0]
        return sequence_data.get(stripped)
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
            filename + ".pickle" for filename in self.fasta_filenames
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

    def _load_or_create_fasta_dictionary_pickle(self):
        self._fasta_dictionary = dict()
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
