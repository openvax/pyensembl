"""
Templates for URLs and paths to specific relase, species, and file type
on the Ensembl ftp server.

For example, the human chromosomal DNA sequences for release 77 are in:

    ftp://ftp.ensembl.org/pub/release-78/fasta/homo_sapiens/dna/

"""
from os.path import join

from species import normalize_species_name
from release_info import which_human_reference, check_release_version

ENSEMBL_FTP_SERVER = 'ftp://ftp.ensembl.org'

# Example directories
# FASTA files: /pub/release-78/fasta/homo_sapiens/
# GTF annotation files: /pub/release-78/gtf/homo_sapiens/
SPECIES_SUBDIR_TEMPLATE = '/pub/release-%(release)d/%(filetype)s/%(species)s/'

# GTF annotation file example: Homo_sapiens.GTCh38.gtf.gz
GTF_FILENAME_TEMPLATE = "%(Species)s.%(reference)s.%(release)d.gtf.gz"

# DNA fasta file example: Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
FASTA_FILENAME_TEMPLATE = \
    "%(Species)s.%(reference)s.%(sequence_type)s.chromosome.%(contig)s.fa.gz"

def _construct_url_subdir(ensembl_release, species):
    """
    Assume ensembl_release has already been normalize by calling function
    but species might be either a common name or latin name.
    """
    subdir = SPECIES_SUBDIR_TEMPLATE % {
        'release' : ensembl_release,
        'filetype' : 'gtf',
        'species' : species,
    }
    return join(ENSEMBL_FTP_SERVER, subdir)

def _normalize_release_properties(ensembl_release, species):
    """
    Make sure a given release is valid, normalize it to be an integer,
    normalize the species name, and get its associated reference.
    """
    ensembl_release = check_release_version(ensembl_release)
    species = normalize_species_name(species)
    # TODO: generalize this to species other than human
    assert species == 'homo_sapiens'
    reference_name = which_human_reference(ensembl_release)
    return ensembl_release, species, reference_name

def gtf_url(ensembl_release, species):
    ensembl_release, species, reference_name = _normalize_release_properties(
        ensembl_release, species)
    url_subdir = _construct_url_subdir(ensembl_release, species)

    filename = GTF_FILENAME_TEMPLATE % {
        'Species' : species.capitalize(),
        'reference' : reference,
        'release' : ensembl_release,
    }
    return join(url_subdir, filename)

def fasta_url(ensembl_release, species, contig, sequence_type="dna"):
    ensembl_release, species, reference_name = _normalize_release_properties(
        ensembl_release, species)
    url_subdir = _construct_url_subdir(ensembl_release, species)
    sequence_subdir = join(url_subdir, sequence_type)
    filename = FASTA_FILENAME_TEMPLATE % {
        "Species" : species.capitalize(),
        "reference" : reference_name,
        "sequence_type" : sequence_type,
        "contig" : contig
    }
    return join(sequence_subdir, filename)


    }
