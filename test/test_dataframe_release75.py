from pyensembl import EnsemblRelease

ensembl = EnsemblRelease(75, auto_download=True)

def test_release_75_length():
    df = ensembl.gtf.dataframe()
    assert df is not None
    # TODO(tavi) This length seems to be different
    # depending on when the GTF was downloaded. Investigate
    # more.
    assert len(df) == 2828312 or len(df) == 2828317

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

