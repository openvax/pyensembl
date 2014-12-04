"""
Check that pyensembl returns correct exon objects for exon IDs
and loci. Make sure the information on the exon object matches
the expected gene ID and location.
"""

from pyensembl import EnsemblRelease

from test_common import releases

ensembl77 = EnsemblRelease(77)

def test_exon_object_by_id_ensembl77_CTNNB1():
    exon = ensembl77.exon_by_id("ENSE00003464041")
    assert exon.gene_name == "CTNNB1", exon.transcript_name
    assert exon.contig == "3", exon.contig
    assert exon.strand == "+"
    assert exon.start ==  41224526
    assert exon.end == 41224753
    assert exon.length ==  228


def test_exon_object_by_locus_ensembl77_CTNNB1():
    exons = ensembl77.exons_at_locus(3, 41224526, strand="+")
    for exon in exons:
        assert exon.gene_name == "CTNNB1", exon.transcript_name
        assert exon.contig == "3", exon.contig
        assert exon.strand == "+"
        assert exon.start <=  41224526
        assert exon.end >= 41224526
