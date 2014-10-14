from os.path import join, exists

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

class EnsemblRelease(object):
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

        # lazily cache DataFrame with expanded attribute columns
        self._local_csv_path = None

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

    def local_csv_path(self):
        """
        Path to CSV which the annotation data with expanded columns
        for optional attributes
        """
        if self._local_csv_path is None:
            gtf_path = self.local_gtf_path()
            if gtf_path.endswith(".gtf"):
                base = gtf_path[:-4]
            else:
                assert gtf_path.endswith(".gtf.gz"), \
                    "Invalid file extension for GTF: %s" % gtf_path
                base = gtf_path[:-7]
            self._local_csv_path = base + ".column_attributes.csv"
        return self._local_csv_path

    def _load_from_gtf(self):
        """
        Parse this release's GTF file and load it as a Pandas DataFrame
        """
        path = self.local_gtf_path()
        return load_gtf_as_dataframe(path)

    def _load_from_csv(self):
        """
        If we've already saved the DataFrame as a CSV
        (with the attributes field expanded into columns),
        then load it. Otherwise parse the GTF file, and save it
        as a CSV
        """
        csv_path = self.local_csv_path()
        if exists(csv_path):
            print "Reading Dataframe from %s" % csv_path
            return pd.read_csv(csv_path)
        else:
            df = self._load_from_gtf()
            df.to_csv(csv_path)
            return df

    def dataframe(self):
        if self._df is None:
            self._df = self._load_from_csv()
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




