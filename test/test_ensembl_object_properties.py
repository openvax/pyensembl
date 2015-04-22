"""
Tests for simple properties of an EnsemblRelease object which don't
require database lookups.
"""

from __future__ import absolute_import

from nose.tools import eq_
from pyensembl import EnsemblRelease

def test_reference_name():
    eq_(EnsemblRelease(release=54).reference_name, "NCBI36")
    eq_(EnsemblRelease(release=74).reference_name, "GRCh37")
    eq_(EnsemblRelease(release=75).reference_name, "GRCh37")
    eq_(EnsemblRelease(release=78).reference_name, "GRCh38")
    eq_(EnsemblRelease(release=79).reference_name, "GRCh38")