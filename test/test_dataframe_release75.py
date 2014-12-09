from pyensembl import EnsemblRelease

ensembl = EnsemblRelease(75)

def test_release_75_length():
    df = ensembl.gtf.dataframe()
    assert df is not None
    assert len(df) == 2828312

def test_release_75_contig_MT():
    df = ensembl.gtf.dataframe(contig="M")
    assert df is not None
    assert len(df) > 0
    assert (df.seqname == "MT").all()
    df = ensembl.gtf.dataframe(contig="MT")
    assert df is not None
    assert len(df) > 0
    assert (df.seqname == "MT").all()


def test_release_75_contig_1():
    df = ensembl.gtf.dataframe(contig="1")
    assert df is not None
    assert len(df) > 0
    assert (df.seqname == "1").all()
    df = ensembl.gtf.dataframe(contig=1)
    assert df is not None
    assert len(df) > 0
    assert (df.seqname == "1").all()

