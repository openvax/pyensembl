from pyensembl import EnsemblRelease, MAX_ENSEMBL_RELEASE

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

def test_max_ensembl_release():
    assert isinstance(MAX_ENSEMBL_RELEASE, int), \
        "Unexpected type for MAX_ENSEMBL_RELEASE: %s" % (
            type(MAX_ENSEMBL_RELEASE),)
    assert 83 <= MAX_ENSEMBL_RELEASE < 1000, \
        "Unexpected value for MAX_ENSEMBL_RELEASE: %d" % MAX_ENSEMBL_RELEASE

def test_int_version():
    for version in range(54, MAX_ENSEMBL_RELEASE):
        EnsemblRelease(version)

def test_str_version():
    for version in range(54, MAX_ENSEMBL_RELEASE):
        EnsemblRelease(str(version))
