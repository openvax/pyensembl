"""Make sure we're getting correct transcritp sequence from Ensembl and that
it's a sequence type which correctly implements `complement`
and `reverse_complement`
"""

from __future__ import absolute_import
from .common import eq_
from pyensembl import genome_for_reference_name

grch38 = genome_for_reference_name("GRCh38")


def test_transcript_sequence_ensembl_grch38():
    # extremely short TRD gene
    seq = grch38.transcript_sequence("ENST00000448914")
    expected = "ACTGGGGGATACG"
    eq_(seq, expected)
    # now try via a Transcript object
    eq_(grch38.transcript_by_id("ENST00000448914").sequence, expected)


def test_coding_sequence_matches_position_ranges():
    # Regression test for https://github.com/openvax/pyensembl/issues/176:
    # len(coding_sequence) used to exceed the total length of
    # coding_sequence_position_ranges by exactly 3 because Ensembl's CDS
    # feature excludes the stop codon while coding_sequence includes it.
    for transcript_id in [
        "ENST00000311936",
        "ENST00000371085",
        "ENST00000275493",
    ]:
        transcript = grch38.transcript_by_id(transcript_id)
        cds_length = len(transcript.coding_sequence)
        ranges_length = sum(
            end - start + 1 for (start, end) in transcript.coding_sequence_position_ranges
        )
        eq_(cds_length, ranges_length)
