"""Make sure we're getting correct transcritp sequence from Ensembl and that
it's a sequence type which correctly implements `complement`
and `reverse_complement`
"""

from __future__ import absolute_import
from nose.tools import eq_
from pyensembl import cached_release

ensembl54 = cached_release(54)
ensembl83 = cached_release(83)

def test_transcript_sequence_ensembl54():
    seq = ensembl54.transcript_sequence("ENST00000321606")
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
    eq_(str(seq), full_transcript_sequence)

    # now get the same sequence via a Transcript object
    eq_(ensembl54.transcript_by_id("ENST00000321606").sequence, seq)


def test_transcript_sequence_ensembl83():
    # extremely short TRD gene
    seq = ensembl83.transcript_sequence("ENST00000448914")
    expected = "ACTGGGGGATACG"
    eq_(seq, expected)
    # now try via a Transcript object
    eq_(ensembl83.transcript_by_id("ENST00000448914").sequence, expected)
