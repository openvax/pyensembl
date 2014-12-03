from pyensembl import EnsemblRelease

from test_common import releases

ensembl77 = EnsemblRelease(77)

def test_exon_object_by_id_ensembl77():
    exon = ensembl77.exon_by_id("ENSE00003464041")
    assert exon.transcript_name == "CTNNB1-001", exon.transcript_name
    assert exon.contig == "3", exon.contig
    assert exon.strand == "+"
    assert exon.start ==  41224526
    assert exon.end == 41224753
    assert exon.frame == 1
    assert exon.length ==  228


def test_exon_object_by_locus_ensembl77():
    exons = ensembl77.exons_by_locus(3, 41224526, strand="+")
    assert len(exons) == 1
    exon = exons[0]
    assert exon.transcript_name == "CTNNB1-001", exon.transcript_name
    assert exon.contig == "3", exon.contig
    assert exon.strand == "+"
    assert exon.start ==  41224526
    assert exon.end == 41224753
    assert exon.frame == 1
    assert exon.length ==  228
