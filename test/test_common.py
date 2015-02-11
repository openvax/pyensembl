import functools

from pyensembl import EnsemblRelease
from nose.tools import nottest

_release_versions = [
    54, # last release for GRCh36/hg18
    75, # last release for GRCh37/hg19
    78, # most recent release for GRCh38
]

_cached_releases = {}

def cached_release(version):
    if version in _cached_releases:
        return _cached_releases[version]
    ensembl = EnsemblRelease(version, auto_download=True)
    _cached_releases[version] = ensembl
    return ensembl

releases = [
    cached_release(version) for version in _release_versions
]

contigs = list(range(1,23)) + ["X", "Y", "M"]

@nottest
def test_ensembl_releases(*versions):
    """
    Run a unit test which takes an EnsemblRelease as an argument
    for multiple rleleases (most recent for each reference genome)
    """
    if len(versions) == 0:
        versions = _release_versions
    def decorator(test_fn):
        @functools.wraps(test_fn)
        def new_test_fn():
            for version in versions:
                ensembl = cached_release(version)
                test_fn(ensembl)
        return new_test_fn
    return decorator
