from pyensembl import EnsemblRelease
from test_common import test_ensembl_releases

#mapping of ensembl releases to transcript IDs for FOXP3-001
FOXP3_transcript_id = "ENST00000376207"

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
    transcript = ensembl.transcript_by_id(FOXP3_transcript_id)
    # The combined lengths of the upstream untranslated region,
    # coding sequence, and downstream untranslated region
    full_sequence = transcript.sequence
    utr5 = transcript.five_prime_utr_sequence
    cds = transcript.coding_sequence
    utr3 = transcript.three_prime_utr_sequence
    # need to use `seq` property of Sequence objects to get underlying
    # strings which can be concatenated and compared
    combined_string = utr5.seq + cds.seq + utr3.seq
    assert combined_string == full_sequence.seq, \
        "Expected FOXP3-001 sequence:\n%s\n\n5' UTR + CDS + 3' UTR:\n%s" % (
            full_sequence,
            combined_string
        )

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

SMAD4_transcript_id = "ENST00000398417"

if __name__ == "__main__":
    test_sequence_parts()
