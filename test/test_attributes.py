
from pyensembl import EnsemblRelease

def test_all_entries_have_gene_id():
    data = EnsemblRelease(75)
    df = data.dataframe()
    assert 'gene_id' in df
    # Ensembl gene ids are formatted like ENSG00000223972
    # which is always length 15
    assert (df['gene_id'].str.len() == 15).all(), \
        df[df['gene_id'].str.len() != 15]
