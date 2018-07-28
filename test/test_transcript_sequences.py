"""Make sure we're getting correct transcritp sequence from Ensembl and that
it's a sequence type which correctly implements `complement`
and `reverse_complement`
"""

from __future__ import absolute_import
from nose.tools import eq_
from pyensembl import genome_for_reference_name

grch38 = genome_for_reference_name("GRCh38")

def test_transcript_sequence_ensembl_grch38():
    # extremely short TRD gene
    seq = grch38.transcript_sequence("ENST00000448914")
    expected = "ACTGGGGGATACG"
    eq_(seq, expected)
    # now try via a Transcript object
    eq_(grch38.transcript_by_id("ENST00000448914").sequence, expected)
