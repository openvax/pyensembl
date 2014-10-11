from os.path import join

from gtf import load_gtf_as_dataframe
from locus import normalize_chromosome

import datacache
import pandas as pd

MIN_ENSEMBL_RELEASE = 48
MAX_ENSEMBL_RELEASE = 77



def _check_release(release):
    """
    Convert a user-provided release number into
    an integer, check to make sure it's in the
    valid range of Ensembl releases
    """
    try:
        release = int(release)
    except:
       assert False, "%s is not a valid Ensembl release" % release
    assert release >= MIN_ENSEMBL_RELEASE
    assert release <= MAX_ENSEMBL_RELEASE
    return release



# mapping from Ensembl release to which reference assembly it uses
_human_references = {}

# Ensembl release 48-54 use NCBI36 as a reference
for i in xrange(48,55):
    _human_references[i] = 'NCBI36'

# Ensembl releases 55-75 use CRCh37 as a reference
for i in xrange(55,76):
    _human_references[i] = 'GRCh37'

# Ensembl releases 76 and 77 use GRCh38
for i in xrange(76,78):
    _human_references[i] = 'GRCh38'

def _which_human_reference(release):
    release = _check_release(release)
    assert release in _human_references, \
        "No reference found for release %d" % release
    return _human_references[release]


# directory which contains GTF files, missing the release number
URL_DIR_TEMPLATE = 'ftp://ftp.ensembl.org/pub/release-%d/gtf/homo_sapiens/'
FILENAME_TEMPLATE = "Homo_sapiens.%s.%d.gtf.gz"

class HumanData(object):
    def __init__(self, release):
        self.release = _check_release(release)
        self.gtf_url_dir = URL_DIR_TEMPLATE % self.release
        self.reference_name =  _which_human_reference(self.release)
        self.gtf_filename = FILENAME_TEMPLATE  % (
            self.reference_name, self.release
        )
        self.gtf_url = join(self.gtf_url_dir, self.gtf_filename)

        # lazily download GTF data if anything beyond path/URL is necessary
        self._local_gtf_path = None

        # lazily load DataFrame if necessary
        self._df = None


    def local_gtf_path(self):
        """
        Returns local path to GTF file for given release of Ensembl,
        download from the Ensembl FTP server if not already cached.
        """
        if self._local_gtf_path is None:
            self._local_gtf_path = datacache.fetch_file(
                self.gtf_url,
                filename=self.gtf_filename,
                decompress=False,
                subdir="ensembl")
        assert self._local_gtf_path
        return self._local_gtf_path

    def dataframe(self):
        if self._df is None:
            path = self.local_gtf_path()
            df = load_gtf_as_dataframe(path)
        assert self._df is not None
        return self._df

    def genes_at_locus(self, chromosome, position):
        df = self.dataframe()
        chromosome = normalize_chromosome(chromosome)
        df_genes = df[df.feature == 'gene']
        df_chr = df_genes[df_genes.seqname == chromosome]

        # find genes whose start/end boundaries overlap with the position
        overlap_start = df_chr.start <= position
        overlap_end = df_chr.end >= position
        overlap = overlap_start & overlap_end
        df_overlap = df_chr[overlap]

        # making a set to only keep unique genes
        genes = set([])
        # parse gene_id from semi-colon separated attribute list

        return genes




