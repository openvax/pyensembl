from __future__ import absolute_import

import functools

from pyensembl import EnsemblRelease
from nose.tools import nottest


_cached_releases = {}

def cached_release(version):
    if version in _cached_releases:
        return _cached_releases[version]
    ensembl = EnsemblRelease(version, auto_download=True)
    _cached_releases[version] = ensembl
    return ensembl

ensembl_grch36 = cached_release(54)  # last release for GRCh36/hg18
ensembl_grch37 = cached_release(75)  # last release for GRCh37/hg19
ensembl_grch38 = cached_release(78)  # most recent release for GRCh38

major_releases = [
    ensembl_grch37,
    ensembl_grch38
]

contigs = list(range(1, 23)) + ["X", "Y", "M"]

@nottest
def test_ensembl_releases(*versions):
    """
    Run a unit test which takes an EnsemblRelease as an argument
    for multiple releases (most recent for each reference genome)
    """
    if len(versions) == 0:
        ensembl_releases = major_releases
    else:
        ensembl_releases = [cached_release(version) for version in versions]

    def decorator(test_fn):
        @functools.wraps(test_fn)
        def new_test_fn():
            for ensembl in ensembl_releases:
                test_fn(ensembl)
        return new_test_fn
    return decorator
