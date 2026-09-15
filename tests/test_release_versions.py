import pytest

from pyensembl import EnsemblRelease, MAX_ENSEMBL_RELEASE
from pyensembl.ensembl_versions import check_release_number

from pytest import raises


@pytest.mark.parametrize("release", [(), (81,), (81, 82), "invalid", None])
def test_invalid_release_preserves_validation_error(release):
    with raises(ValueError) as error:
        check_release_number(release)
    assert str(error.value) == "Invalid Ensembl release: %s" % (release,)


@pytest.mark.parametrize("release", [81, "81"])
def test_valid_release_normalization(release):
    assert check_release_number(release) == 81


def test_version_too_old_1():
    with raises(Exception):
        EnsemblRelease(1)


def test_version_too_old_47():
    with raises(Exception):
        EnsemblRelease(47)


def test_version_is_not_numeric():
    with raises(Exception):
        EnsemblRelease("wuzzle")


def test_version_is_none():
    with raises(Exception):
        EnsemblRelease(None)


def test_max_ensembl_release():
    assert isinstance(
        MAX_ENSEMBL_RELEASE, int
    ), "Unexpected type for MAX_ENSEMBL_RELEASE: %s" % (type(MAX_ENSEMBL_RELEASE),)
    assert 83 <= MAX_ENSEMBL_RELEASE < 1000, (
        "Unexpected value for MAX_ENSEMBL_RELEASE: %d" % MAX_ENSEMBL_RELEASE
    )


def test_int_version():
    for version in range(54, MAX_ENSEMBL_RELEASE):
        EnsemblRelease(version)


def test_str_version():
    for version in range(54, MAX_ENSEMBL_RELEASE):
        EnsemblRelease(str(version))
