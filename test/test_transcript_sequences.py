"""Make sure we're getting correct transcritp sequence from Ensembl and that
it's a sequence type which correctly implements `complement`
and `reverse_complement`
"""

from __future__ import absolute_import
from nose.tools import eq_
from pyensembl import ensembl_grch38

def test_transcript_sequence_ensembl_grch38():
    # extremely short TRD gene
    seq = ensembl_grch38.transcript_sequence("ENST00000448914")
    expected = "ACTGGGGGATACG"
    eq_(seq, expected)
    # now try via a Transcript object
    eq_(ensembl_grch38.transcript_by_id("ENST00000448914").sequence, expected)
