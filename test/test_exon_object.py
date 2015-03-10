"""
Check that pyensembl returns correct exon objects for exon IDs
and loci. Make sure the information on the exon object matches
the expected gene ID and location.
"""
from __future__ import absolute_import

from pyensembl import EnsemblRelease

ensembl77 = EnsemblRelease(77, auto_download=True)

def test_exon_object_by_id():
    """
    test_exon_object_by_id : check properties of exon 4 of CTNNB1 when looked
    up by ID in Ensembl 77.
    """
    exon = ensembl77.exon_by_id("ENSE00003464041")
    assert exon.gene_name == "CTNNB1", \
        "Unexpected gene name: %s" % exon.gene_name
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
    assert exon.gene_name == "CXCR3", \
        "Unexpected gene name: %s" % exon.gene_name
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
    transcript = ensembl77.transcripts_by_name("CXCR3-002")[0]
    exon = transcript.exons[0]
    assert exon.contains_start_codon

def test_contains_stop_codon():
    """
    Test that second exon of CXCR3-002 contains a stop codon.
    """
    transcript = ensembl77.transcripts_by_name("CXCR3-002")[0]
    exon = transcript.exons[1]
    assert exon.contains_stop_codon

def test_partial_start_codons():
    """
    test_partial_start_codons : Make sure none of the start codons returned
    only partially overlap with an exon.
    """
    # Exon ENSE00003718948 has a start codon immediately before the
    # exon's start position on the forward strand
    exon = ensembl77.exon_by_id("ENSE00003718948")
    for (start, end) in exon.start_codon_positions:
        assert start >= exon.start, \
            "Exon at locus [%d, %d], start codon at %d" % (
                exon.start, exon.end, start)

    for offset in exon.start_codon_offsets:
        assert offset >= 0, \
            "Invalid negative offset for start codon: %d" % offset
        assert offset != 0, \
            "Partially overlapping start codon should not be included"

def test_start_codon_offset():
    """
    test_start_codon_offset : Ensure that start codon in exon #1 of CXCR3-002
    is 68 nucleotides from the start of the exon.
    """
    transcript = ensembl77.transcripts_by_name("CXCR3-002")[0]
    exon = transcript.exons[0]
    start_offsets = exon.start_codon_offsets
    assert len(start_offsets) == 1, \
        "First exon of CXCR3-002 should only overlap one start codon"
    offset = start_offsets[0]
    assert offset == 68, \
        "Expected stop codon 68nt from start of exon #1, got %s" % (offset,)

def test_stop_codon_offset():
    """
    test_stop_codon_offset : Ensure that the only stop codon
    overlapping exon #2 of CXCR-002 is 1092 nucleotides
    from the start of the exon.
    """
    transcript = ensembl77.transcripts_by_name("CXCR3-002")[0]
    exon = transcript.exons[1]
    stop_offsets = exon.stop_codon_offsets
    assert len(stop_offsets) == 1, \
        "Last exon of CXCR3-002 should only overlap one stop codon, got: %s" % (
        stop_offsets)
    offset = stop_offsets[0]
    assert offset == 1092, \
        "Expected stop codon 1092nt from start of exon #2, got %s" % (offset,)

