"""
The datacache module already gives us download caching,
this goes a step further by also caching parsed CSVs in-memory
and serializing the results of expensive computations to CSVs
"""

from os.path import exists

import pandas as pd

def _read_csv(csv_path):
    print "Reading Dataframe from %s" % csv_path
    df = pd.read_csv(csv_path)
    # by default, Pandas will infer the type as int,
    # then switch to str when it hits non-numerical
    # chromosomes. Make sure whole column has the same type
    df['seqname'] = df['seqname'].map(str)
    return df

def _write_csv(df, csv_path):
    print "Saving DataFrame to %s" % csv_path
    df.to_csv(csv_path, index=False)

_dataframe_cache = {}

def clear_cached_objects():
    if len(_dataframe_cache) > 0:
        print "Removing %d items from DataFrame cache" % len(_dataframe_cache)
        _dataframe_cache.clear()

def remove_from_cache(key):
    if key in _dataframe_cache:
        del _dataframe_cache[key]

def load_csv(csv_path, expensive_action):
    """
    If we've already saved the DataFrame as a CSV
    (with the attributes field expanded into columns),
    then load it. Otherwise parse the GTF file, and save it
    as a CSV
    """
    if csv_path in _dataframe_cache:
        return _dataframe_cache[csv_path]

    if exists(csv_path):
        df = _read_csv(csv_path)
    else:
        df = expensive_action()
        _write_csv(df, csv_path)
    _dataframe_cache[csv_path] = df
    return df