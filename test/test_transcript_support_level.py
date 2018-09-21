"""
Tests for methods which return collections of transcript IDs that aren't
converting from some type of name or ID.
"""
from __future__ import absolute_import

from nose.tools import eq_

from pyensembl import cached_release

def test_transcript_support_level():
    ensembl93 = cached_release(93)
    transcript = ensembl93.transcripts_by_name("DDX11L1-202")[0]
    eq_(transcript.support_level, "1")
    # for this transcript, the tsl value is missing, and should yield an empty string:
    transcript = ensembl93.transcripts_by_name("OR4G11P-202")[0]
    eq_(transcript.support_level, '')

    # transcript_support_level columns was missing in GRCh37 and older releases of GRCh38, and should yield None:
    ensembl77 = cached_release(77)
    transcript = ensembl77.transcripts_by_name("DDX11L1-002")[0]
    eq_(transcript.support_level, None)
