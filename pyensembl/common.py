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

from functools import wraps

CACHE_SUBDIR = "ensembl"

def memoize(fn):
    """Simple memoization decorator for functions and methods,
    assumes that all arguments to the function can be hashed and
    compared.
    """
    memoized_values = {}

    @wraps(fn)
    def wrapped_fn(*args, **kwargs):
        key = (args, tuple(sorted(kwargs.items())))
        if key not in memoized_values:
            memoized_values[key] = fn(*args, **kwargs)
        return memoized_values[key]

    return wrapped_fn
