from os.path import join

from common import CACHE_SUBDIR
from release_info import which_human_reference_name, check_release_number
from url_templates import ENSEMBL_FTP_SERVER, fasta_cdna_url_parts

import pyfaidx
from datacache import Cache

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
            server=ENSEMBL_FTP_SERVER):

        # download cache for fetching reference FASTA files
        self.cache = Cache(CACHE_SUBDIR)

        self.release = check_release_number(ensembl_release)

        assert species == "homo_sapiens", \
            "Species other than human not supported: %s" % (species,)
        self.reference_name = which_human_reference_name(ensembl_release)
        assert species == "homo_sapiens", \
            "Species %s not currently supported" % (species,)
        self.species = species

        self.server = server

        reference_url_dir, reference_filename = fasta_cdna_url_parts(
            ensembl_release=self.release,
            species=self.species,
            server=server)
        self.remote_filename = reference_filename
        self.url = join(reference_url_dir, reference_filename)

        # dictionary mapping transcript IDs to cDNA sequences
        self._transcript_sequences = {}

        self._fasta_dictionary = None

    @property
    def local_fasta_path(self):
        """
        Returns local path to FASTA file,
        download from the Ensembl FTP server if not already cached.
        """
        return self.cache.fetch(self.url, self.remote_filename, decompress=True)

    @property
    def local_dir(self):
        return split(self.local_fasta_path)[0]

    @property
    def local_filename(self):
        return split(self.local_fasta_path)[1]

    @property
    def fasta_dictionary(self):
        if not self._fasta_dictionary:
            self._fasta_dictionary = pyfaidx.Fasta(self.local_fasta_path)
        return self._fasta_dictionary

    def transcript_sequence(self, transcript_id):
        if transcript_id not in self._transcript_sequences:
            if not transcript_id.startswith("ENST"):
                raise ValueError("Invalid transcript ID: %s" % (transcript_id,))

            if transcript_id not in self.fasta_dictionary:
                raise ValueError(
                    "Transcript ID not found: %s" % (transcript_id,))
            fasta_record = self.fasta_dictionary[transcript_id]

            # FastaRecord doesn't seem to have an accessor to get the full
            # sequence (only subsequences), so slice out the full string
            seq = fasta_record[:len(fasta_record)]
            self._transcript_sequences[transcript_id] = seq
        return self._transcript_sequences[transcript_id]

    def __str__(self):
        return "ReferenceTranscripts(release=%s, species=%s, filename=%s)"  % (
            self.release, self.species, self.remote_filename)

    def __repr__(self):
        return str(self)

    def __contains__(self, transcript_id):
        return transcript_id in self.fasta_dictionary
