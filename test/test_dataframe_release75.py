from __future__ import absolute_import

from pyensembl import EnsemblRelease

ensembl = EnsemblRelease(75, auto_download=True)


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

