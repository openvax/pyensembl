from ensembl.human_data import EnsemblRelease

def test_release_75_length():
    ensembl = EnsemblRelease(75)
    df = ensembl.dataframe()
    assert df is not None
    assert len(df) == 2828312

def test_release_75_contig_MT():
    df = ensembl.dataframe_for_contig("M")
    assert df
    assert len(df) > 0
    assert (df.seqname == "MT").all()
    df = ensembl.dataframe_for_contig("MT")
    assert df
    assert len(df) > 0
    assert (df.seqname == "MT").all()


def test_release_75_contig_1():
    df = ensembl.dataframe_for_contig("1")
    assert df
    assert len(df) > 0
    assert (df.seqname == "MT").all()
    df = ensembl.dataframe_for_contig(1)
    assert df
    assert len(df) > 0
    assert (df.seqname == "MT").all()

