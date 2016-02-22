"""
Tests for simple properties of an EnsemblRelease object which don't
require database lookups.
"""

from __future__ import absolute_import

from nose.tools import eq_
from pyensembl import EnsemblRelease, MAX_ENSEMBL_RELEASE

def test_human_reference_name():
    eq_(EnsemblRelease(release=54).reference_name, "NCBI36")
    eq_(EnsemblRelease(release=74).reference_name, "GRCh37")
    eq_(EnsemblRelease(release=75).reference_name, "GRCh37")
    for release in range(76, MAX_ENSEMBL_RELEASE):
        eq_(EnsemblRelease(release=release).reference_name, "GRCh38")
