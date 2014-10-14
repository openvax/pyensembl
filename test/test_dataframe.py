from ensembl.human_data import EnsemblRelease

def test_release_75_length():
    ensembl = EnsemblRelease(75)
    df = ensembl.dataframe()
    assert df is not None
    assert len(df) == 2828312
