from os.path import exists

from ensembl.human_data import EnsemblRelease

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
    data = EnsemblRelease(75)
    path = data.local_gtf_path()
    assert exists(path)
    assert path.endswith(".gtf.gz")

def test_version_75_string():
    data = EnsemblRelease("75")
    path = data.local_gtf_path()
    assert exists(path)
    assert path.endswith(".gtf.gz")
