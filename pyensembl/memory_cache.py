# Copyright (c) 2015. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Cache and serializing the results of expensive computations. Used in pyensembl
primarily to cache the heavy-weight parsing of GTF files and various
filtering operations on Ensembl entries.

A piece of data is returned from one of three sources:
1) Cache cold. Run the user-supplied compute_fn.
2) Cache warm on disk. Parse or unpickle the serialized result into memory.
3) Cache warm in memory. Return cached object.
"""
from __future__ import print_function, division, absolute_import

import logging
from os import remove, stat
from os.path import exists

import pandas as pd

from .common import load_pickle, dump_pickle

class MemoryCache(object):
    def __init__(self):
        self._memory_cache = {}

    def is_empty(self, filename):
        return stat(filename).st_size == 0

    def delete_file(self, path):
        if exists(path):
            logging.info("Deleting cached file %s" % path)
            remove(path)

    def remove_from_cache(self, key):
        if key in self._memory_cache:
            del self._memory_cache[key]
        self.delete_file(key)

    def clear_cached_objects(self):
        for key in self._memory_cache.keys():
            self.delete_file(key)
        self._memory_cache.clear()

    def _read_csv(self, csv_path):
        print("Reading Dataframe from %s" % csv_path)
        df = pd.read_csv(csv_path)
        if 'seqname' in df:
            # by default, Pandas will infer the type as int,
            # then switch to str when it hits non-numerical
            # chromosomes. Make sure whole column has the same type
            df['seqname'] = df['seqname'].map(str)
        return df

    def _write_csv(self, df, csv_path, chunksize=10**5):
        """
        Parameters
        ----------
        df : pandas.DataFrame

        csv_path : str

        chunksize : int
            Number of rows to write at a time. Helps to limit memory
            consumption while writing a CSV.
        """
        print("Saving DataFrame to %s" % csv_path)
        df.to_csv(csv_path, index=False, chunksize=chunksize)

    def cached_dataframe(self, csv_path, compute_fn):
        """
        If a CSV path is in the _memory_cache, then return that cached value.

        If we've already saved the DataFrame as a CSV then load it.

        Otherwise run the provided `compute_fn`, and store its result
        in memory and and save it as a CSV.
        """
        if not csv_path.endswith(".csv"):
            raise ValueError("Invalid path '%s', must be a CSV file" % csv_path)

        if csv_path in self._memory_cache:
            return self._memory_cache[csv_path]

        if exists(csv_path) and not self.is_empty(csv_path):
            df = self._read_csv(csv_path)
        else:
            df = compute_fn()
            if not isinstance(df, pd.DataFrame):
                raise TypeError(
                    "Expected compute_fn to return DataFrame, got %s : %s" % (
                        df, type(df)))
            self._write_csv(df, csv_path)
        self._memory_cache[csv_path] = df
        return df

    def cached_object(self, path, compute_fn):
        """
        If `cached_object` has already been called for a value of `path` in this
        running Python instance, then it should have a cached value in the
         _memory_cache; return that value.

        If this function was never called before with a particular value of
        `path`, then call compute_fn, and pickle it to `path`.

        If `path` already exists, unpickle it and store that value in
        _memory_cache.
        """
        if path in self._memory_cache:
            return self._memory_cache[path]

        if exists(path) and not self.is_empty(path):
            obj = load_pickle(path)
        else:
            obj = compute_fn()
            dump_pickle(obj, path)
        self._memory_cache[path] = obj
        return obj
