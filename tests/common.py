import functools

from pyensembl import genome_for_reference_name, cached_release

import pytest


grch37 = genome_for_reference_name("GRCh37")
grch38 = genome_for_reference_name("GRCh38")

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


# TemporaryDirectory only got added to Python in version 3.2
try:
    # pylint: disable=no-name-in-module
    from tempfile import TemporaryDirectory

except ImportError:
    # only added in Python 3.2
    from tempfile import mkdtemp
    from shutil import rmtree

    class TemporaryDirectory(object):
        def __init__(self):
            self.name = mkdtemp()

        def __enter__(self, *args, **kwargs):
            return self.name

        def __exit__(self, type, value, traceback):
            rmtree(self.name)
            # don't suppress exceptions
            return False


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
