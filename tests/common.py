
from pyensembl import cached_release

import pytest

# The releases CI installs (.github/workflows/tests.yml), pinned so other
# releases in the local cache cannot change what the tests check (#398).
grch37 = cached_release(75)
grch38 = cached_release(93)

major_releases = [grch37, grch38]

contigs = [str(c) for c in range(1, 23)] + ["X", "Y", "M"]


def run_multiple_genomes(*versions):
    if len(versions) == 1 and callable(versions[0]):
        return pytest.mark.parametrize("genome", major_releases)(versions[0])
    if not versions:
        genomes = major_releases
    else:
        genomes = [cached_release(v) for v in versions]
    return lambda fn: pytest.mark.parametrize("genome", genomes)(fn)


def ok_(b):
    assert b


def eq_(x, y, msg=None):
    if msg is None:
        assert x == y
    else:
        assert x == y, msg


def neq_(x, y, msg=None):
    if msg is None:
        assert x != y
    else:
        assert x != y, msg


def gt_(x, y, msg=None):
    if msg is None:
        assert x > y
    else:
        assert x > y, msg


def gte_(x, y, msg=None):
    if msg is None:
        assert x >= y
    else:
        assert x >= y, msg


def lt_(x, y, msg=None):
    if msg is None:
        assert x < y
    else:
        assert x < y, msg


def lte_(x, y, msg=None):
    if msg is None:
        assert x <= y
    else:
        assert x <= y, msg
