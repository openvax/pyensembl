"""
Test the ReferenceTranscripts object to make sure we're getting
the correct reference for different Ensembl releases and that the
the correct transcripts get returned for known transcript IDs.
"""

from pyensembl import EnsemblRelease

def test_reference_name():
    reference54 = EnsemblRelease(release=54).reference
    assert reference54.reference_name == "NCBI36"

    reference75 = EnsemblRelease(release=75).reference
    assert reference75.reference_name == "GRCh37"

    reference77 = EnsemblRelease(release=77).reference
    assert reference77.reference_name == "GRCh38"

def test_transcript_sequence_ensembl54():
    reference54 = EnsemblRelease(release=54,
                                 auto_download=True).reference
    seq = reference54.transcript_sequence("ENST00000321606")
    assert len(seq) == 414, \
        "Expected transcript ENST00000321606 to have 414nt, got %s : %d" % (
            seq, len(seq))
    nucleotide_lines = [
        "CATGTCACCCACCTTCAGGCGGCCCAAGACACTGCGACTCCGGAGGCAGCCCAGATATCCTCGGAAGAG",
        "CACCCCCAGGAGAAACAAGCTTGGCCACTATGCTATCATCAAGTTTCCGCTGACCACTGAGTCGGCCGT",
        "GAAGAAGATAGAAGAAAACAACACGCTTGTGTTCACTGTGGATGTTAAAGCCAACAAGCACCAGATCAG",
        "ACAGGCTGTGAAGAAGCTCTATGACAGTGATGTGGCCAAGGTCACCACCCTGATTTGTCCTGATAAAGA",
        "GAACAAGGCATATGTTCGACTTGCTCCTGATTATGATGCTTTCGATGTTGTAACAAAATTGGGATCACC",
        "TAAACTGAGTCCAGCTGGCTAACTCTAAATATATGTGTATCTTTTCAGCATAAAAAAATAATGTTTTTC"
    ]
    full_transcript_sequence = "".join(nucleotide_lines)
    assert str(seq) == full_transcript_sequence

    # make sure the Sequence object's reverse and complement properties work
    reverse_complement = seq.complement
    reverse_complement_lines = [
        "GAAAAACATTATTTTTTTATGCTGAAAAGATACACATATATTTAGAGTTAGCCAGCTGGACTCAGTTTA",
        "GGTGATCCCAATTTTGTTACAACATCGAAAGCATCATAATCAGGAGCAAGTCGAACATATGCCTTGTTC",
        "TCTTTATCAGGACAAATCAGGGTGGTGACCTTGGCCACATCACTGTCATAGAGCTTCTTCACAGCCTGT",
        "CTGATCTGGTGCTTGTTGGCTTTAACATCCACAGTGAACACAAGCGTGTTGTTTTCTTCTATCTTCTTC",
        "ACGGCCGACTCAGTGGTCAGCGGAAACTTGATGATAGCATAGTGGCCAAGCTTGTTTCTCCTGGGGGTG",
        "CTCTTCCGAGGATATCTGGGCTGCCTCCGGAGTCGCAGTGTCTTGGGCCGCCTGAAGGTGGGTGACATG"
    ]
    reverse_complement = "".join(reverse_complement_lines)
    assert len(seq.reverse.complement) == 414
    assert str(seq.reverse.complement) == reverse_complement
