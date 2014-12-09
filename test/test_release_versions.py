from os.path import exists

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

def test_version_75():
    ensembl = EnsemblRelease(75)
    path = ensembl.gtf.local_gtf_path()
    assert exists(path)
    assert path.endswith(".gtf.gz")

def test_version_75_string():
    ensembl = EnsemblRelease("75")
    path = ensembl.gtf.local_gtf_path()
    assert exists(path)
    assert path.endswith(".gtf.gz")
