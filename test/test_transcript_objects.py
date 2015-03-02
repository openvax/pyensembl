from pyensembl import Locus, Transcript

from test_common import cached_release, test_ensembl_releases
from data import (
    FOXP3_001_transcript_id,
    CTNNBIP1_004_transcript_id,
    CTNNBIP1_004_UTR5,
    CTNNBIP1_004_UTR3,
    CTNNBIP1_004_CDS,
    CTNNBIP1_004_locus,
    CTTNNIP1_004_exon_lengths,
    CTTNNIP1_004_exon_ids,
)

def test_transcript_start_codon():
    """
    test_transcript_start_codon : Check that fields Transcript
    (for transcript named CTNNBIP1-004) matches known values.
    """
    ensembl = cached_release(77)
    CTNNBIP1_004_transcript = ensembl.transcript_by_id(
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
    ensembl = cached_release(77)
    transcript = ensembl.transcript_by_id(CTNNBIP1_004_transcript_id)
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
                i+1, transcript, expected_id, exon.id)

        expected_length = CTTNNIP1_004_exon_lengths[i]
        assert len(exon) == expected_length, \
            "Expected exon #%d of %s (%s) to have length %d but got %d" % (
                i+1, transcript, exon, expected_length, len(exon))


# not testing NCBI/Release 54 since I just discovered that ensembl54
# feature='transcript' entries don't have a gene ID.
# TODO: Add gene_id patching to gtf_parsing, add ensembl54 to the list
# below
@test_ensembl_releases(75, 77)
def test_sequence_parts(ensembl):
    """
    test_sequence_parts : Ensure that the UTRs and coding sequence can be
    combined to make the full transcript.
    """
    transcript = ensembl.transcript_by_id(FOXP3_001_transcript_id)

    # The combined lengths of the upstream untranslated region,
    # coding sequence, and downstream untranslated region
    full_sequence = transcript.sequence
    assert len(full_sequence) > 0

    utr5 = transcript.five_prime_utr_sequence
    assert len(utr5) > 0

    cds = transcript.coding_sequence
    assert len(cds) > 0

    utr3 = transcript.three_prime_utr_sequence
    assert len(utr3) > 0

    # need to use `seq` property of Sequence objects to get underlying
    # strings which can be concatenated and compared
    combined_string = utr5.seq + cds.seq + utr3.seq

    combined_sequence_length = len(combined_string)
    # make sure length property of transcript matches the sequence length
    assert combined_sequence_length == len(transcript), \
        "Length 5' UTR(%dnt) + CDS(%dnt) + 3' UTR(%d) = %d, expected %d" % (
            len(utr5),
            len(cds),
            len(utr3),
            combined_sequence_length,
            len(transcript)
        )
    assert combined_string == full_sequence.seq, \
        "Expected FOXP3-001 sequence:\n%s\n\n5' UTR + CDS + 3' UTR:\n%s" % (
            full_sequence,
            combined_string
        )

def test_transcript_sequences_CTNNIP1_004():
    ensembl = cached_release(77)
    transcript = ensembl.transcript_by_id(CTNNBIP1_004_transcript_id)
    utr5 = transcript.five_prime_utr_sequence
    expected_utr5_length = len(CTNNBIP1_004_UTR5)
    assert len(utr5) == expected_utr5_length, \
        "Expected 5' UTR length %d, got %d" % (
            expected_utr5_length, len(utr5))
    assert utr5.seq == CTNNBIP1_004_UTR5, "5' UTR sequence is incorrect"

    utr3 = transcript.three_prime_utr_sequence
    expected_utr3_length = len(CTNNBIP1_004_UTR3)
    assert len(utr3) == expected_utr3_length, \
        "Expected 3' UTR length %d, got %d" % (
            expected_utr3_length, len(utr3))
    assert utr3.seq == CTNNBIP1_004_UTR3, "3' UTR sequence is incorrect"

    cds = transcript.coding_sequence
    expected_cds_length = len(CTNNBIP1_004_CDS)
    assert len(cds) == expected_cds_length, \
        "Expected CDS length %d, got %d" % (expected_cds_length, len(cds))
    assert cds.seq == CTNNBIP1_004_CDS, "Coding sequence is incorrect"


@test_ensembl_releases
def test_equal_transcripts(release):
    t1 = release.transcripts_by_name("TP53-001")[0]
    # make an identical gene
    t2 = Transcript(t1.id, t1.db, t1.reference)

    assert hash(t1) == hash(t2)
    assert t1 == t2

@test_ensembl_releases
def test_not_equal_transcripts(release):
    t1 = release.genes_by_name("MUC1-001")[0]
    t2 = release.genes_by_name("BRCA1-001")[0]
    assert hash(t1) != hash(t2)
    assert t1 != t2
