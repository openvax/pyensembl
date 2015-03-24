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
from os.path import join, exists, split
import logging

from Bio import SeqIO
import datacache

from .common import CACHE_SUBDIR
from .release_info import which_human_reference_name, check_release_number
from .url_templates import ENSEMBL_FTP_SERVER, fasta_cdna_url_parts


class ReferenceTranscripts(object):
    """
    Container for reference genome sequenes. Downloads and caches reference
    FASTA files.
    Currently only supports (unspliced) transcript sequences.

    TODO: Access introns and intergenic regions by also downloading
    the *.dna.fa.gz FASTA files. Should rename class to Reference when
    this gets implemented.
    """
    def __init__(
            self,
            ensembl_release,
            species="homo_sapiens",
            server=ENSEMBL_FTP_SERVER,
            auto_download=False):

        # download cache for fetching reference FASTA files
        self.cache = datacache.Cache(CACHE_SUBDIR)
        self.auto_download = auto_download

        self.release = check_release_number(ensembl_release)

        assert species == "homo_sapiens", \
            "Species other than human not supported: %s" % (species,)
        self.reference_name = which_human_reference_name(ensembl_release)
        assert species == "homo_sapiens", \
            "Species %s not currently supported" % (species,)
        self.species = species

        self.server = server

        self.fasta_decompress = True
        reference_url_dir, reference_filename = fasta_cdna_url_parts(
            ensembl_release=self.release,
            species=self.species,
            server=server)
        self.remote_filename = reference_filename
        self.url = join(reference_url_dir, reference_filename)

        # dictionary mapping transcript IDs to cDNA sequences
        self._transcript_sequences = {}

        # SeqIO.index_db gets lazily constructed since we only want to download
        # the Fasta file if/when we use it
        self._fasta_dictionary = None

        # key set for fasta dictionary, for faster membership tests,
        self._fasta_keys = None

    def __str__(self):
        return "ReferenceTranscripts(release=%s, species=%s, filename=%s)" % (
            self.release, self.species, self.remote_filename)

    def __repr__(self):
        return str(self)

    def __contains__(self, transcript_id):
        if self._fasta_keys is None:
            self._fasta_keys = set(self.fasta_dictionary.keys())
        return transcript_id in self._fasta_keys

    def __eq__(self, other):
        return (
            isinstance(other, ReferenceTranscripts) and
            self.release == other.release and
            self.species == other.species and
            self.server == other.server)

    @property
    def local_fasta_path(self):
        """
        Returns local path to FASTA file. If it's not already
        cached, download it from the Ensembl FTP server if auto
        download is enabled.
        """
        # If the fasta is already cached, fetching it won't initiate a
        # download. But it's always okay to initiate a download if
        # auto download is enabled.
        if (self.cache.exists(
                self.url,
                self.remote_filename,
                self.fasta_decompress) or self.auto_download):
            # Does a download if the cache is empty.
            return self.cache.fetch(self.url,
                                    self.remote_filename,
                                    self.fasta_decompress)
        raise ValueError("Ensembl transcript data is not currently "
                         "installed for release %s. Run "
                         "\"pyensembl install %s\" or call "
                         "EnsemblRelease(%s).install()" %
                         ((self.release,) * 3))

    @property
    def local_dir(self):
        return split(self.local_fasta_path)[0]

    @property
    def local_filename(self):
        return split(self.local_fasta_path)[1]

    def transcript_sequence(self, transcript_id, return_none_if_missing=True):
        if transcript_id not in self._transcript_sequences:
            if not transcript_id.startswith("ENST"):
                raise ValueError("Invalid transcript ID: %s" % (transcript_id,))

            seq_record = self.fasta_dictionary.get(transcript_id)
            if seq_record is None:
                if return_none_if_missing:
                    return None
                raise ValueError(
                    "Transcript ID not found: %s" % (transcript_id,))
            self._transcript_sequences[transcript_id] = seq_record.seq
        return self._transcript_sequences[transcript_id]

    def download_transcript_sequences(self, force=False):
        """
        Download the FASTA file if one does not exist. If `force` is
        True, overwrites any existing file.

        Returns True if a download happened.
        """
        if not force and self.cache.exists(self.url,
                                           self.remote_filename,
                                           self.fasta_decompress):
            return False
        self.cache.fetch(self.url, self.remote_filename,
                         self.fasta_decompress, force=force)
        return True

    @property
    def local_database_path(self):
        return self.local_fasta_path + ".db"

    def _create_or_open_fasta_db(self):
        """Create an index database from a transcript sequence FASTA file"""
        return SeqIO.index_db(
            self.local_database_path,
            self.local_fasta_path,
            "fasta")

    def index(self, force=True):
        """
        Create a database mapping transcript IDs to sequences.

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
                    "Deleting existing transcript database: %s",
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
            # backed dictionary of reference transcript sequences
            self.index()
        return self._fasta_dictionary
