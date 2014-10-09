from os.path import exists

from ensembl.human_data import HumanData

from nose.tools import raises


@raises(Exception)
def test_version_too_old_1():
    HumanData(1)

@raises(Exception)
def test_version_too_old_47():
    HumanData(47)

@raises(Exception)
def test_version_is_not_numeric():
    HumanData("wuzzle")

@raises(Exception)
def test_version_is_none():
    HumanData(None)

def test_version_75():
    data = HumanData(75)
    path = data.local_gtf_path()
    assert exists(path)
    assert path.endswith(".gtf.gz")

def test_version_75_string():
    data = HumanData("75")
    path = data.local_gtf_path()
    assert exists(path)
    assert path.endswith(".gtf.gz")