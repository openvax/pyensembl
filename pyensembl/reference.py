from os.path import join

from common import CACHE_SUBDIR
from release_info import which_human_reference_name
from url_templates import ENSEMBL_FTP_SERVER, fasta_cdna_url

from datacache import Cache

class ReferenceTranscripts(object):
    def __init__(
            self,
            ensembl_release,
            species="homo_sapiens",
            server=ENSEMBL_FTP_SERVER):
        self.cache = Cache(CACHE_SUBDIR)
        self.release = ensembl_release
        assert species == "homo_sapiens", \
            "Species %s not currently supported" % (species,)
        self.species = species
        self.server = server

        reference_url_dir, reference_filename = fasta_cdna_url(
            release=self.release,
            species=self.species,
            server=server)
        self.filename = reference_filename
        self.reference_url = join(reference_url_dir, reference_filename)
