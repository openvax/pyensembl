# Copyright (c) 2015. Mount Sinai School of Medicine
#
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

from __future__ import print_function, division, absolute_import
from os import remove
from os.path import exists, split, isfile
import logging

from Bio import SeqIO
import datacache

from .common import CACHE_SUBDIR, require_ensembl_id


class SequenceData(object):
    """
    Container for reference nucleotide and amino acid sequenes.

    Downloads and caches FASTA file, creates database mapping
    identifiers to sequences.
    """
    def __init__(
            self,
            genome_source,
            fasta_type,
            local_filename_func=None,
            require_ensembl_ids=True,
            auto_download=False):

        # download cache for fetching reference FASTA files
        self.cache = datacache.Cache(CACHE_SUBDIR)
        self.auto_download = auto_download
        # expect the remote FASTA file to be gzipped and the local needs
        # to be a decompressed text file
        self.fasta_decompress = True
        self.path = genome_source.fasta_path(fasta_type)
        self.remote_filename = split(self.path)[1]

        if local_filename_func:
            self.local_filename = local_filename_func(self.remote_filename)
        else:
            self.local_filename = self.remote_filename
        self._init_lazy_fields()

    def _init_lazy_fields(self):
        """Initialize all properties which get lazily constructed on request"""
        # dictionary mapping IDs to sequences
        self._sequence_cache = {}

        # SeqIO.index_db gets lazily constructed since we only want to download
        # the Fasta file if/when we use it
        self._fasta_dictionary = None

        # key set for fasta dictionary, for faster membership tests,
        self._fasta_keys = None

    def __str__(self):
        return "SequenceData(path=%s, version=%s, local_path=%s)" % (
            self.path,
            self.version,
            self.local_fasta_path)

    def __repr__(self):
        return str(self)

    def __contains__(self, sequence_id):
        if self._fasta_keys is None:
            self._fasta_keys = set(self.fasta_dictionary.keys())
        return sequence_id in self._fasta_keys

    def __eq__(self, other):
        return (
            isinstance(other, SequenceData) and
            self.path == other.path and
            self.version == other.version)

    @property
    def local_fasta_path(self):
        """
        Returns local path to FASTA file. If it's not already
        cached, download it from the Ensembl FTP server if auto
        download is enabled.
        """
        # If the path is a local file as opposed to a URL, then
        # that's our local path.
        if isfile(self.path):
            return self.path

        # If the fasta is already cached, fetching it won't initiate a
        # download. But it's always okay to initiate a download if
        # auto download is enabled.
        if self.local_file_exists() or self.auto_download:
            # Does a download if the cache is empty, otherwise just returns
            # the local path
            return self._fetch(force=False)
        raise ValueError("Genome %s sequence data %s is not currently "
                         "installed for this genome source. Run %s "
                         " or call %s" % (
                             self.fasta_type,
                             self.genome_source.install_string_console(),
                             self.genome_source.install_string_python()))

    def clear_cache(self):
        """Delete the local FASTA file and its associated database,
        and clear any in-memory cached sequence data.
        """
        # reset in-memory cached object
        self._init_lazy_fields()

        if self.local_file_exists():
            remove(self.local_fasta_path)
        if exists(self.local_database_path):
            remove(self.local_fasta_path)

    @property
    def local_dir(self):
        return split(self.local_fasta_path)[0]

    def local_file_exists(self):
        return self.cache.exists(
                self.path,
                filename=self.local_filename,
                decompress=self.fasta_decompress)

    def _fetch(self, force):
        """Download file from remote URL if not present or if force=True

        Return local path
        """
        return self.cache.fetch(
            url=self.path,
            filename=self.local_filename,
            decompress=self.fasta_decompress,
            force=force)

    def download(self, force=False):
        """Download the FASTA file if one does not exist.

        If `force` is True, overwrites any existing file.
        """
        # If the path is a local file as opposed to a URL, then
        # no download is necesary.
        if isfile(self.path):
            return

        if not self.local_file_exists() or force:
            self._fetch(force=force)

    @property
    def local_database_path(self):
        return self.local_fasta_path + ".db"

    def _create_or_open_fasta_db(self):
        """Create an index database from a transcript sequence FASTA file"""
        return SeqIO.index_db(
            self.local_database_path,
            self.local_fasta_path,
            "fasta")

    def index(self, force=False):
        """Create a database mapping transcript IDs to sequences.

        Parameters
        -----------
        force : bool
            Recreate database even if it already exists

        Raises an error if the necessary data is not yet downloaded.

        If no error is raised, then after call self._fasta_dictionary
        should be populated.
        """
        # This local_fasta_path property access will raise an error
        # if the necessary data is not yet downloaded
        if exists(self.local_database_path):
            if force:
                logging.info(
                    "Deleting existing sequence database: %s",
                    self.local_database_path)
                remove(self.local_database_path)
            else:
                # try opening the existing FASTA database, if it fails then
                # delete the local file and we'll create a new one further
                # below
                try:
                    self._fasta_dictionary = self._create_or_open_fasta_db()
                    return
                except ValueError as e:
                    # if there was an error opening the database
                    # delete the version we have and try again
                    if "database" in e.message:
                        logging.warn(
                            ("Recreating transcript database %s"
                             " due to error: %s"),
                            self.local_database_path,
                            e.message)
                        remove(self.local_database_path)
                    else:
                        raise
        # if database didn't exist or was incomplete, create it now
        self._fasta_dictionary = self._create_or_open_fasta_db()

    @property
    def fasta_dictionary(self):
        if self._fasta_dictionary is None:
            # should populate self._fasta_dictionary with a database
            # backed dictionary of reference sequences
            self.index()
        return self._fasta_dictionary

    def get(self, sequence_id):
        """Get sequence associated with given ID or return None if missing"""

        if sequence_id not in self._sequence_cache:
            # all Ensembl identifiers start with ENS e.g. ENST, ENSP, ENSE
            if require_ensembl_ids:
                require_ensembl_id(sequence_id)

            # get sequence from database
            fasta_record = self.fasta_dictionary.get(sequence_id)
            if fasta_record:
                seq = fasta_record.seq
            else:
                seq = None
            self._sequence_cache[sequence_id] = seq
        return self._sequence_cache[sequence_id]
