"""
Check that pyensembl returns correct exon objects for exon IDs
and loci. Make sure the information on the exon object matches
the expected gene ID and location.
"""
from __future__ import absolute_import

from pyensembl import ensembl_grch38 as ensembl

def test_exon_object_by_id():
    """
    test_exon_object_by_id : check properties of exon 4 of CTNNB1 when looked
    up by ID in Ensembl 77.
    """
    exon = ensembl.exon_by_id("ENSE00003464041")
    assert exon.gene_name == "CTNNB1", \
        "Unexpected gene name: %s" % exon.gene_name
    assert exon.contig == "3", exon.contig
    assert exon.strand == "+"
    assert exon.on_forward_strand
    assert exon.on_positive_strand
    assert exon.start == 41224526, "Unexpected exon start: %s" % exon.start
    assert exon.end == 41224753, "Unexpected exon end: %s" % exon.end
    assert exon.length == len(exon) == 228

def test_exon_object_by_id_on_negative_strand():
    """
    test_exon_object_by_id : check properties of exon 1 from CXCR3 when looked
    up by ID in Ensembl 77.
    """
    exon = ensembl.exon_by_id("ENSE00001817013")
    assert exon.gene_name == "CXCR3", \
        "Unexpected gene name: %s" % exon.gene_name
    assert exon.contig == "X", exon.contig
    assert exon.strand == "-"
    assert exon.on_backward_strand
    assert exon.on_negative_strand
    assert exon.start == 71618438, "Unexpected exon start: %s" % exon.start
    assert exon.end == 71618517, "Unexpected exon end: %s" % exon.end
    assert exon.length == len(exon) == 80


def test_exon_object_at_locus():
    """
    test_exon_object_at_locus : check properties of exon 4 of CTNNB1 when looked
    up by its location on the forward strand of chr3
    """
    exons = ensembl.exons_at_locus(3, 41224526, strand="+")
    for exon in exons:
        assert exon.gene_name == "CTNNB1", exon.transcript_name
        assert exon.contig == "3", exon.contig
        assert exon.strand == "+"
        assert exon.on_forward_strand
        assert exon.on_positive_strand
        assert exon.start <= 41224526, "Unexpected exon start: %s" % exon.start
        assert exon.end >= 41224526, "Unexpected exon end: %s" % exon.end

def test_exon_object_at_locus_on_negative_strand():
    """
    test_exon_object_at_locus : check properties of exon 1 of CXCR3 when looked
    up by its location on the negative strand of chrX
    """
    exons = ensembl.exons_at_locus("chrX", 71618517, strand="-")
    for exon in exons:
        assert exon.gene_name == "CXCR3", exon.transcript_name
        assert exon.contig == "X", exon.contig
        assert exon.strand == "-"
        assert exon.on_backward_strand
        assert exon.on_negative_strand
        assert exon.start <= 71618517, "Unexpected exon start: %s" % exon.start
        assert exon.end >= 71618517, "Unexpected exon end: %s" % exon.end
