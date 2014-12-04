"""
Check that pyensembl returns correct exon objects for exon IDs
and loci. Make sure the information on the exon object matches
the expected gene ID and location.
"""

from pyensembl import EnsemblRelease

from test_common import releases

ensembl77 = EnsemblRelease(77)

def test_exon_object_by_id():
    """
    test_exon_object_by_id : check properties of exon 4 of CTNNB1 when looked
    up by ID in Ensembl 77.
    """
    exon = ensembl77.exon_by_id("ENSE00003464041")
    assert exon.gene_name == "CTNNB1", exon.transcript_name
    assert exon.contig == "3", exon.contig
    assert exon.strand == "+"
    assert exon.on_forward_strand
    assert exon.on_positive_strand
    assert exon.start ==  41224526, "Unexpected exon start: %s" % exon.start
    assert exon.end == 41224753, "Unexpected exon end: %s" % exon.end
    assert exon.length ==  228

def test_exon_object_by_id_on_negative_strand():
    """
    test_exon_object_by_id : check properties of exon 1 from CXCR3 when looked
    up by ID in Ensembl 77.
    """
    exon = ensembl77.exon_by_id("ENSE00001817013")
    assert exon.gene_name == "CXCR3", exon.transcript_name
    assert exon.contig == "X", exon.contig
    assert exon.strand == "-"
    assert exon.on_backward_strand
    assert exon.on_negative_strand
    assert exon.start ==  71618438, "Unexpected exon start: %s" % exon.start
    assert exon.end == 71618517, "Unexpected exon end: %s" % exon.end
    assert exon.length ==  80


def test_exon_object_at_locus():
    """
    test_exon_object_at_locus : check properties of exon 4 of CTNNB1 when looked
    up by its location on the forward strand of chr3
    """
    exons = ensembl77.exons_at_locus(3, 41224526, strand="+")
    for exon in exons:
        assert exon.gene_name == "CTNNB1", exon.transcript_name
        assert exon.contig == "3", exon.contig
        assert exon.strand == "+"
        assert exon.on_forward_strand
        assert exon.on_positive_strand
        assert exon.start <=  41224526, "Unexpected exon start: %s" % exon.start
        assert exon.end >= 41224526, "Unexpected exon end: %s" % exon.end

def test_exon_object_at_locus_on_negative_strand():
    """
    test_exon_object_at_locus : check properties of exon 1 of CXCR3 when looked
    up by its location on the negative strand of chrX
    """
    exons = ensembl77.exons_at_locus("chrX", 71618517, strand="-")
    for exon in exons:
        assert exon.gene_name == "CXCR3", exon.transcript_name
        assert exon.contig == "X", exon.contig
        assert exon.strand == "-"
        assert exon.on_backward_strand
        assert exon.on_negative_strand
        assert exon.start <=  71618517, "Unexpected exon start: %s" % exon.start
        assert exon.end >= 71618517, "Unexpected exon end: %s" % exon.end


def test_contains_start_codon():
    """
    test_contains_start_codon : Test that first exon of CXCR3-002 contains
    a start codon.
    """
    transcript = ensembl77.transcript_by_name("CXCR3-002")
    print transcript
    print transcript.exons
    exon = transcript.exons[0]
    assert exon.contains_start_codon

def test_contains_stop_codon():
    """
    Test that second exon of CXCR3-002 contains a stop codon.
    """
    transcript = ensembl77.transcript_by_name("CXCR3-002")
    exon = transcript.exons[1]
    assert exon.contains_stop_codon

def test_start_codon_offset():
    """
    test_start_codon_offset : Ensure that start codon in exon #1 of CXCR3-002
    is 68 nucleotides from the start of the exon.
    """
    transcript = ensembl77.transcript_by_name("CXCR3-002")
    exon = transcript.exons[0]
    offset = exon.start_codon_offset
    assert offset == 68, \
        "Expected stop codon 68nt from start of exon #1, got %s" % (offset,)

def test_stop_codon_offset():
    """
    test_stop_codon_offset : Ensure that stop codon in exon #2 of CXCR-002
    is 1098 nucleotides from the start of the exon.
    """
    transcript = ensembl77.transcript_by_name("CXCR3-002")
    print transcript
    print transcript.exons
    exon = transcript.exons[1]
    offset = exon.stop_codon_offset
    assert offset == 1092, \
        "Expected stop codon 1092nt from start of exon #2, got %s" % (offset,)

