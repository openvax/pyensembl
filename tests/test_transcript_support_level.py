"""
Tests for methods which return collections of transcript IDs that aren't
converting from some type of name or ID.
"""
from __future__ import absolute_import

from .common import eq_

from pyensembl import cached_release
from pyensembl.genome import _parse_transcript_support_level


def test_parse_transcript_support_level_values():
    # Recent Ensembl releases append parenthetical text to the TSL; we should
    # still recover the leading integer. See openvax/pyensembl#297.
    eq_(_parse_transcript_support_level("1 (assigned to previous version 5)"), 1)
    eq_(_parse_transcript_support_level("3 (assigned to previous version 12)"), 3)
    # Plain numeric strings still work.
    eq_(_parse_transcript_support_level("1"), 1)
    eq_(_parse_transcript_support_level("5"), 5)
    # Missing / NA / junk values collapse to None.
    eq_(_parse_transcript_support_level(None), None)
    eq_(_parse_transcript_support_level(""), None)
    eq_(_parse_transcript_support_level("NA"), None)
    eq_(_parse_transcript_support_level("NA (something)"), None)
    eq_(_parse_transcript_support_level("   "), None)


def test_transcript_support_level():
    """The Transcript Support Level (TSL) is a method to highlight the well-supported and poorly-supported transcript
    models for users, based on the type and quality of the alignments used to annotate the transcript.
    In the Ensembl database, it can be assigned to a value 1 through 5, or reported as NA, or missing, or missing
    completely in older releases. We translate it to an integer value, otherwise to None.
    """
    ensembl93 = cached_release(93)
    transcript = ensembl93.transcripts_by_name("DDX11L1-202")[0]
    eq_(transcript.support_level, 1)

    # For this transcript, the transcript_support_level value is missing in the database record:
    transcript = ensembl93.transcripts_by_name("OR4G11P-202")[0]
    eq_(transcript.support_level, None)

    # Some features are reported as "NA" in Ensembl: those are features like pseudogenes, single exon transcripts,
    # HLA, T-cell receptor and Ig transcripts that are not analysed in terms of TSL and therefore not given any
    # of the TSL categories. We translate NA to None as well.
    transcript = ensembl93.transcripts_by_name("MIR1302-2-201")[0]
    eq_(transcript.support_level, None)

    # Transcript_support_level column was missing completely in GRCh37 and older releases of GRCh38:
    ensembl77 = cached_release(77)
    transcript = ensembl77.transcripts_by_name("DDX11L1-002")[0]
    eq_(transcript.support_level, None)
