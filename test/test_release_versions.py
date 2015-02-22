from pyensembl import EnsemblRelease

from nose.tools import raises


@raises(Exception)
def test_version_too_old_1():
    EnsemblRelease(1)

@raises(Exception)
def test_version_too_old_47():
    EnsemblRelease(47)

@raises(Exception)
def test_version_is_not_numeric():
    EnsemblRelease("wuzzle")

@raises(Exception)
def test_version_is_none():
    EnsemblRelease(None)

def test_int_version():
    for version in range(54, 77):
        EnsemblRelease(version)

def test_str_version():
    for version in range(54, 77):
        EnsemblRelease(str(version))
