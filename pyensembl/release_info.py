from __future__ import print_function, division, absolute_import

MIN_ENSEMBL_RELEASE = 48
MAX_ENSEMBL_RELEASE = 78

def check_release_number(release):
    """
    Convert a user-provided release number into
    an integer, check to make sure it's in the
    valid range of Ensembl releases
    """
    try:
        release = int(release)
    except:
        raise ValueError("Invalid Ensembl release: %s" % release)

    if release < MIN_ENSEMBL_RELEASE or release > MAX_ENSEMBL_RELEASE:
        raise ValueError(
            "Invalid Ensembl releases %d, must be between %d and %d" %
                (release, MIN_ENSEMBL_RELEASE, MAX_ENSEMBL_RELEASE))
    return release

# mapping from Ensembl release to which reference assembly it uses
_human_references = {}

# Ensembl release 48-54 use NCBI36 as a reference
for i in range(48,54+1):
    _human_references[i] = 'NCBI36'

# Ensembl releases 55-75 use GRCh37 as a reference
for i in range(55,75+1):
    _human_references[i] = 'GRCh37'

# Most recent Ensembl releases use GRCh38
for i in range(76,MAX_ENSEMBL_RELEASE+1):
    _human_references[i] = 'GRCh38'

def which_human_reference_name(release):
    release = check_release_number(release)
    if release not in _human_references:
        raise ValueError(
            "No reference found for release %d" % release)
    return _human_references[release]
