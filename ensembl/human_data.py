import datacache
import pandas as pd

MIN_ENSEMBL_RELEASE = 48
MAX_ENSEMBL_RELEASE = 77

def check_release(release):
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


# directory which contains GTF files, missing the release number
URL_DIR_TEMPLATE = 'ftp://ftp.ensembl.org/pub/release-%d/gtf/homo_sapiens/'

def gtf_url_dir(release):
    release = check_release(release)
    return URL_DIR_TEMPLATE % release


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

def which_human_reference(release):
    release = check_release(release)
    assert release in _human_references, \
        "No reference found for release %d" % release
    return _human_references[release]

FILENAME_TEMPLATE = "Homo_sapiens.%s.%d.gtf.gz"

def gtf_filename(release):
    release = check_release(release)
    reference_name = which_human_reference(release)
    return FILENAME_TEMPLATE % (reference_name, release)


def gtf_url(release):
    base_url = gtf_url_dir(release)
    filename = gtf_filename(release)
    return base_url + filename

def local_gtf_path(release):
    """
    Returns local path to GTF file for given release of Ensembl,
    download from the Ensembl FTP server if not already cached.
    """
    filename = gtf_filename(release)
    url = gtf_url(release)
    local_path = datacache.fetch_file(
        url,
        filename=filename,
        decompress=True,
        subdir="ensembl")
    return local_path

GTF_COLS = [
    'seqname',
    'source',
    'feature',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'attribute'
]

def load_dataframe(release):
    """
    Download GTF annotation data for given release (if not already present),
    parse it into a dataframe.
    """
    path = local_gtf_path(release)
    return pd.read_csv(path, comment='#', sep='\t', names = GTF_COLS)

