from __future__ import absolute_import

from pyensembl import ensembl_grch37


def test_release_75_contig_MT():
    df = ensembl_grch37.gtf.dataframe(contig="M")
    assert df is not None
    assert len(df) > 0
    assert (df.seqname == "MT").all()
    df = ensembl_grch37.gtf.dataframe(contig="MT")
    assert df is not None
    assert len(df) > 0
    assert (df.seqname == "MT").all()


def test_release_75_contig_1():
    df = ensembl_grch37.gtf.dataframe(contig="1")
    assert df is not None
    assert len(df) > 0
    assert (df.seqname == "1").all()
    df = ensembl_grch37.gtf.dataframe(contig=1)
    assert df is not None
    assert len(df) > 0
    assert (df.seqname == "1").all()
