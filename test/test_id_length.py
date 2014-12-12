from test_common import releases

from nose.tools import nottest
from pyensembl import EnsemblRelease

@nottest
def check_id_length(feature_name):
    for release in releases:
        # only load chromosome Y to speed up tests
        df = release.gtf.dataframe(contig="Y")
        assert feature_name in df, \
            "%s not a column in DataFrame for %s" % (feature_name, release)
        # Ensembl IDs are formatted like ENSG00000223972
        # which is always length 15
        ids = df[feature_name]
        # any entry missing this ID should be dropped
        id_lengths = ids.str.len().dropna()
        correct_mask = id_lengths == 15
        assert correct_mask.all(), \
           "Invalid %s IDs: %s" % (
                feature_name,
                df[~correct_mask][feature_name])

def test_gene_id_length():
    check_id_length('gene_id')

def test_transcript_id_length():
    check_id_length('transcript_id')

def test_protein_id_length():
    check_id_length('protein_id')




