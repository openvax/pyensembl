from __future__ import absolute_import

from pyensembl import Locus, cached_release
from nose.tools import eq_, assert_not_equal, assert_greater

from .common import test_ensembl_releases
from .data import (
    FOXP3_001_transcript_id,
    CTNNBIP1_004_transcript_id,
    CTNNBIP1_004_UTR5,
    CTNNBIP1_004_UTR3,
    CTNNBIP1_004_CDS,
    CTNNBIP1_004_locus,
    CTTNNIP1_004_exon_lengths,
    CTTNNIP1_004_exon_ids,
    EGFR_001_protein_sequence,
    TP53_gene_id,
)

ensembl77 = cached_release(77)

def test_transcript_start_codon():
    """
    test_transcript_start_codon : Check that fields Transcript
    (for transcript named CTNNBIP1-004) matches known values.
    """
    CTNNBIP1_004_transcript = ensembl77.transcript_by_id(
        CTNNBIP1_004_transcript_id)
    assert Locus.__eq__(CTNNBIP1_004_locus, CTNNBIP1_004_transcript), \
        "Expected locus %s but got %s" % (
            CTNNBIP1_004_locus, Locus.__str__(CTNNBIP1_004_transcript))

    start_offsets = CTNNBIP1_004_transcript.start_codon_spliced_offsets
    assert len(start_offsets) == 3, \
        "Wrong length for start codon: %d (%s)" % (
            len(start_offsets), start_offsets)

    assert all(isinstance(i, int) for i in start_offsets), \
        "Wrong type %s for beginning start codon offset" % (
            [type(i) for i in start_offsets],)

    expected_start_codon_offset = len(CTNNBIP1_004_UTR5)
    start_codon_offset = min(start_offsets)
    assert start_codon_offset == expected_start_codon_offset, \
        "Incorrect start codon offset, expected %d but got %d" % (
            expected_start_codon_offset, start_codon_offset)

def test_transcript_exons():
    """
    test_transcript_exons : Ensure that properties of CTTNBIP1-004's
    Exon objects match known values.
    """
    transcript = ensembl77.transcript_by_id(CTNNBIP1_004_transcript_id)
    exons = transcript.exons
    assert isinstance(exons, list), \
        "Expected list of Exon objects, got %s : %s" % (exons, type(exons))

    # CTTNBIP1-004 has 5 exons
    assert len(exons) == len(CTTNNIP1_004_exon_lengths), \
        "Expected %d exons but got %d" % (
            len(CTTNNIP1_004_exon_lengths), len(exons))

    for i, exon in enumerate(exons):
        expected_id = CTTNNIP1_004_exon_ids[i]
        assert exon.id == expected_id, \
            "Expected exon #%d of %s to have ID %s but got %s" % (
                i + 1, transcript, expected_id, exon.id)

        expected_length = CTTNNIP1_004_exon_lengths[i]
        assert len(exon) == expected_length, \
            "Expected exon #%d of %s (%s) to have length %d but got %d" % (
                i + 1, transcript, exon, expected_length, len(exon))


# not testing NCBI/Release 54 since I just discovered that ensembl54
# feature='transcript' entries don't have a gene ID.
# TODO: Add gene_id patching to gtf_parsing, add ensembl54 to the list
# below
@test_ensembl_releases(75, 77)
def test_sequence_parts(genome):
    # Ensure that the UTRs and coding sequence can be
    # combined to make the full transcript.
    transcript = genome.transcript_by_id(FOXP3_001_transcript_id)

    # The combined lengths of the upstream untranslated region,
    # coding sequence, and downstream untranslated region
    full_sequence = transcript.sequence
    assert_greater(len(full_sequence), 0)

    utr5 = transcript.five_prime_utr_sequence
    assert_greater(len(utr5), 0)

    cds = transcript.coding_sequence
    assert_greater(len(cds), 0)

    utr3 = transcript.three_prime_utr_sequence
    assert_greater(len(utr3), 0)

    # need to use `seq` property of Sequence objects to get underlying
    # strings which can be concatenated and compared
    combined_string = utr5 + cds + utr3

    combined_sequence_length = len(combined_string)
    # make sure length property of transcript matches the sequence length
    eq_(
        combined_sequence_length,
        len(transcript),
        "Length 5' UTR(%dnt) + CDS(%dnt) + 3' UTR(%d) = %d, expected %d" % (
            len(utr5),
            len(cds),
            len(utr3),
            combined_sequence_length,
            len(transcript)))
    eq_(
        combined_string,
        full_sequence,
        "Expected FOXP3-001 sequence:\n%s\n\n5' UTR + CDS + 3' UTR:\n%s" % (
            full_sequence,
            combined_string))

def test_transcript_utr5_sequence_CTNNIP1_004():
    transcript = ensembl77.transcript_by_id(CTNNBIP1_004_transcript_id)
    utr5 = transcript.five_prime_utr_sequence
    expected_utr5_length = len(CTNNBIP1_004_UTR5)
    eq_(len(utr5),
        expected_utr5_length,
        "Expected 5' UTR length %d, got %d" % (
            expected_utr5_length, len(utr5)))
    eq_(utr5, CTNNBIP1_004_UTR5)

def test_transcript_utr3_sequence_CTNNIP1_004():
    transcript = ensembl77.transcript_by_id(CTNNBIP1_004_transcript_id)
    utr3 = transcript.three_prime_utr_sequence
    expected_utr3_length = len(CTNNBIP1_004_UTR3)
    eq_(len(utr3),
        expected_utr3_length,
        "Expected 3' UTR length %d, got %d" % (
            expected_utr3_length, len(utr3)))
    eq_(utr3, CTNNBIP1_004_UTR3)

def test_transcript_cds_CTNNIP1_004():
    transcript = ensembl77.transcript_by_id(CTNNBIP1_004_transcript_id)
    cds = transcript.coding_sequence
    expected_cds_length = len(CTNNBIP1_004_CDS)
    eq_(
        len(cds),
        expected_cds_length,
        "Expected CDS length %d, got %d" % (expected_cds_length, len(cds)))
    eq_(cds, CTNNBIP1_004_CDS)

@test_ensembl_releases()
def test_equal_transcripts(genome):
    t1 = genome.transcripts_by_name("TP53-001")[0]
    # get an identical gene
    t2 = genome.transcript_by_id(t1.id)
    eq_(t1, t2)
    eq_(hash(t1), hash(t2))

@test_ensembl_releases()
def test_not_equal_transcripts(genome):
    t1 = genome.transcripts_by_name("MUC1-001")[0]
    t2 = genome.transcripts_by_name("BRCA1-001")[0]
    assert_not_equal(t1, t2)

def test_protein_id():
    transcript = ensembl77.transcripts_by_name("EGFR-001")[0]
    eq_(transcript.protein_id, "ENSP00000275493")

def test_protein_protein_sequence():
    transcript = ensembl77.transcripts_by_name("EGFR-001")[0]
    eq_(transcript.protein_sequence, EGFR_001_protein_sequence)

def test_transcript_gene_should_match_parent_gene():
    gene = ensembl77.gene_by_id(TP53_gene_id)
    for transcript in gene.transcripts:
        eq_(transcript.gene, gene)

@test_ensembl_releases()
def test_BRCA1_001_has_protein_coding_biotype(genome):
    transcript = genome.transcripts_by_name("BRCA1-001")[0]
    assert transcript.is_protein_coding
    eq_(transcript.biotype, "protein_coding")
