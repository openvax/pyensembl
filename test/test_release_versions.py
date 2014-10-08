from os.path import exists

from ensembl.human_data import local_gtf_path

from nose.tools import raises


@raises(Exception)
def test_version_too_old_1():
    local_gtf_path(1)

@raises(Exception)
def test_version_too_old_47():
    local_gtf_path(47)

@raises(Exception)
def test_version_is_not_numeric():
    local_gtf_path("wuzzle")

@raises(Exception)
def test_version_is_none():
    local_gtf_path(None)

def test_version_75():
    path = local_gtf_path(75)
    assert exists(path)
    assert path.endswith(".gtf")

def test_version_75_string():
    path = local_gtf_path("75")
    assert exists(path)
    assert path.endswith(".gtf")