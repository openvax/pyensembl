"""Make sure we're getting correct transcritp sequence from Ensembl and that
it's a sequence type which correctly implements `complement`
and `reverse_complement`
"""

from __future__ import absolute_import
from nose.tools import eq_
from pyensembl import ensembl54

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
