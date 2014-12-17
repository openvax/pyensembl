from os import remove
from os.path import exists
import tempfile

from pyensembl import compute_cache

class Counter(object):
    """
    Use this class to count how many times a function gets called by
    cached_object and cached_dataframe.
    """
    def __init__(self):
        self.count = 0

    def increment(self):
        self.count += 1
        return self.count


def test_cached_object_empty_tempfile():
    """
    test_cached_dataframe_file_exists : A temporary file exists before
    calling into compute_cache but is empty (and thus shouldn't be loaded
    as a null result.
    """
    counter = Counter()
    with tempfile.NamedTemporaryFile() as f:
        first_result = compute_cache.cached_object(
            f.name, compute_fn=counter.increment)
        assert first_result == 1, "Expected count=1, got %s" % (first_result,)
        assert counter.count == 1, \
            "compute_cache.cached_object never called its compute_fn"
        second_result = compute_cache.cached_object(
            f.name, compute_fn=counter.increment)
        assert second_result == 1, \
            "Expected cached count=1, got %s" % (second_result,)
        assert counter.count == 1, \
            "cached_object should have called compute_fn once, got %s" % (
                counter.count)

