from __future__ import absolute_import

import tempfile

from pyensembl import MemoryCache

import pandas as pd
from nose.tools import raises

memory_cache = MemoryCache()

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

    def increment_dataframe(self):
        value = self.increment()
        return pd.DataFrame({'x': [value]})

def test_cached_object_with_tempfile():
    """
    test_cached_object_with_tempfile : A temporary file exists before
    calling into compute_cache.cached_object but is empty, should be treated
    as if result has never been computed before (rather than trying to load
    the empty file).
    """
    counter = Counter()
    with tempfile.NamedTemporaryFile() as f:
        # call repeatedly to test the hot and cold cache logic
        result = memory_cache.cached_object(
            f.name, compute_fn=counter.increment)
        assert result == 1, "Expected result=1, got %s" % (result,)
        assert counter.count == 1, \
            "Expected compute_fn to be called once, got %s" % (counter.count,)


def test_cached_dataframe_with_tempfile():
    """
    test_cached_dataframe_with_tempfile : A temporary file exists before
    calling into compute_cache.cached_dataframe but is empty,
    should be treated as if result has never been computed before
    (rather than trying to load the empty file).
    """
    counter = Counter()
    with tempfile.NamedTemporaryFile(suffix='.csv') as f:
        # call repeatedly to test hot and cold cache logic
        for _ in range(2):
            df = memory_cache.cached_dataframe(
                f.name, compute_fn=counter.increment_dataframe)
            # get counter value from inside of dataframe
            result = df['x'].ix[0]
            assert result == 1, \
                "Expected result=1, got %s" % (result,)
            assert counter.count == 1, \
                "Expected compute_fn to be called once, got %s" % (
                    counter.count,)

def test_cached_dataframe_returns_correct_type():
    def make_a_dataframe():
        return pd.DataFrame({'x': [0, 1, 2]})
    with tempfile.NamedTemporaryFile(suffix='.csv') as f:
        # call repeatedly to test the cold and hot cache logic
        for _ in range(2):
            df = memory_cache.cached_dataframe(
                f.name, compute_fn=make_a_dataframe)
            assert isinstance(df, pd.DataFrame), \
                "Expected DataFrame, got %s : %s" % (df, type(df))

def test_cached_object_with_list_returns_correct_type():
    def make_a_list():
        return [1, 2, 3]
    with tempfile.NamedTemporaryFile() as f:
        # call repeatedly to test the cold and hot cache logic
        for _ in range(2):
            df = memory_cache.cached_object(
                f.name, compute_fn=make_a_list)
            assert isinstance(df, list), \
                "Expected list, got %s : %s" % (df, type(df))

@raises(Exception)
def test_dataframe_path_must_be_csv():
    # compute_cache should raise an exception when filename doesn't
    # end with .csv extension
    memory_cache.cached_dataframe(
        csv_path="tempfile_not_csv",
        compute_fn=lambda _: pd.DataFrame({'x': []}))
