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
from os.path import exists, split
import logging

from Bio import SeqIO
import datacache

from .common import (CACHE_SUBDIR, require_ensembl_id)


class SequenceData(object):
    """
    Container for reference nucleotide and amino acid sequenes.

    Downloads and caches FASTA file, creates database mapping
    identifiers to sequences.

    fasta_source is a GenomeSource object that represents a local
    file path or remote URL.
    """
    def __init__(
            self,
            fasta_source,
            require_ensembl_ids=True,
            auto_download=False):
        self.fasta_source = fasta_source
        # download cache for fetching reference FASTA files
        self.cache = datacache.Cache(CACHE_SUBDIR)
        self.auto_download = auto_download
        # expect the original FASTA file to be gzipped and the cached needs
        # to be a decompressed text file
        self.fasta_decompress = True
        assert self.fasta_source, (
            "Attempted to create SequenceData with missing FASTA source")
        self.require_ensembl_ids = require_ensembl_ids
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
        return "SequenceData(fasta_source=%s, cached_path=%s)" % (
            self.fasta_source,
            self.cached_fasta_path)

    def __repr__(self):
        return str(self)

    def __contains__(self, sequence_id):
        if self._fasta_keys is None:
            self._fasta_keys = set(self.fasta_dictionary.keys())
        return sequence_id in self._fasta_keys

    def __eq__(self, other):
        return (
            isinstance(other, SequenceData) and
            self.fasta_source == other.fasta_source)

    def __hash__(self):
        return hash(self.fasta_source)

    @property
    def cached_fasta_path(self):
        """
        Returns cached path to FASTA file. If it's not already
        cached, download it from the appropriate server (or copy
        from local) if auto download is enabled.
        """
        # If the fasta is already cached, fetching it won't initiate a
        # download. But it's always okay to initiate a download if
        # auto download is enabled.
        if self.cached_file_exists() or self.auto_download:
            # If the source is a local file as opposed to a URL, then
            # just manually copy it to the datacache directory if
            # we're auto-downloading.
            if not self.fasta_source.is_url_format():
                self.fasta_source.copy_to_cache_if_needed(self.cache,
                                                          force=False)

            # Does a download if the cache is empty, otherwise just returns
            # the cached path
            return self._fetch(force=False)
        raise ValueError("Genome sequence data (%s) is not currently "
                         "installed for this genome source. Run %s "
                         "or call %s" % (
                             self.fasta_source,
                             self.fasta_source.install_string_console(),
                             self.fasta_source.install_string_python()))

    def clear_cache(self):
        """Delete the cached FASTA file and its associated database,
        and clear any in-memory cached sequence data.
        """
        # reset in-memory cached object
        self._init_lazy_fields()

        if self.cached_file_exists():
            remove(self.cached_fasta_path)
        if exists(self.cached_database_path):
            remove(self.cached_fasta_path)

    @property
    def cached_dir(self):
        return split(self.cached_fasta_path)[0]

    def cached_file_exists(self):
        # If the source is a local file as opposed to a URL, then
        # check to see if it was copied into the cache directory.
        if not self.fasta_source.is_url_format():
            return exists(self.fasta_source.cached_path(self.cache))

        return self.cache.exists(
            self.fasta_source.path_or_url,
            filename=self.fasta_source.cached_filename,
            decompress=self.fasta_decompress)

    def _fetch(self, force):
        """Download file from remote URL if not present or if force=True

        Return cached path
        """
        return self.cache.fetch(
            url=self.fasta_source.path_or_url,
            filename=self.fasta_source.cached_filename,
            decompress=self.fasta_decompress,
            force=force)

    def download(self, force=False):
        """Download the FASTA file if one does not exist.

        If `force` is True, overwrites any existing file.
        """
        # If the source is a local file as opposed to a URL, then
        # just copy it to the datacache direcotry as a "download".
        if not self.fasta_source.is_url_format():
            return self.fasta_source.copy_to_cache_if_needed(
                self.cache, force=force)

        if not self.cached_file_exists() or force:
            self._fetch(force=force)

    @property
    def cached_database_path(self):
        return self.cached_fasta_path + ".db"

    def _create_or_open_fasta_db(self):
        """Create an index database from a transcript sequence FASTA file"""
        return SeqIO.index_db(
            self.cached_database_path,
            self.cached_fasta_path,
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
        # This cached_fasta_path property access will raise an error
        # if the necessary data is not yet downloaded
        if exists(self.cached_database_path):
            if force:
                logging.info(
                    "Deleting existing sequence database: %s",
                    self.cached_database_path)
                remove(self.cached_database_path)
            else:
                # try opening the existing FASTA database, if it fails then
                # delete the cached file and we'll create a new one further
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
                            self.cached_database_path,
                            e.message)
                        remove(self.cached_database_path)
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
            if self.require_ensembl_ids:
                require_ensembl_id(sequence_id)

            # get sequence from database
            fasta_record = self.fasta_dictionary.get(sequence_id)
            if fasta_record:
                seq = fasta_record.seq
            else:
                seq = None
            self._sequence_cache[sequence_id] = seq
        return self._sequence_cache[sequence_id]
