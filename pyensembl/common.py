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

from contextlib import contextmanager
from functools import wraps
import os
from pathlib import Path
import pickle
from uuid import uuid4


def _staging_path(directory, name):
    # Dot files in cache directories are always staging files.
    return Path(directory) / (".%s.%s.tmp" % (name, uuid4().hex))


def _remove(path):
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass


def _publish(staged, path):
    """Durably replace path with a complete staged file."""
    with open(staged, "r+b") as handle:
        os.fsync(handle.fileno())
    os.replace(staged, path)
    try:
        directory = os.open(os.path.dirname(path), os.O_RDONLY)
    except OSError:
        return  # Directories cannot be opened for fsync on some platforms.
    try:
        os.fsync(directory)
    except OSError:
        pass
    finally:
        os.close(directory)


@contextmanager
def _atomic_output(path, mode="wb"):
    """Write a file that appears complete or not at all.

    The staging file is created like any new file (0o666 less umask), so a
    group-shared cache stays readable by other users.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    staged = _staging_path(path.parent, path.name)
    try:
        descriptor = os.open(staged, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666)
        with open(descriptor, mode) as handle:
            yield handle
        _publish(staged, path)
    finally:
        _remove(staged)


def dump_pickle(obj, filepath):
    # Atomic, so an interrupted write never leaves a truncated index behind.
    with _atomic_output(filepath) as f:
        # use lower protocol for compatibility between Python 2 and Python 3
        pickle.dump(obj, file=f, protocol=2)


def load_pickle(filepath):
    with open(filepath, "rb") as f:
        obj = pickle.load(f)
    return obj


def _memoize_cache_key(args, kwargs):
    """Turn args tuple and kwargs dictionary into a hashable key.

    Expects that all arguments to a memoized function are either hashable
    or can be uniquely identified from type(arg) and repr(arg).
    """
    cache_key_list = []

    # hack to get around the unhashability of lists,
    # add a special case to convert them to tuples
    for arg in args:
        if type(arg) is list:
            cache_key_list.append(tuple(arg))
        else:
            cache_key_list.append(arg)
    for k, v in sorted(kwargs.items()):
        if type(v) is list:
            cache_key_list.append((k, tuple(v)))
        else:
            cache_key_list.append((k, v))
    return tuple(cache_key_list)


def merge_intervals(ranges):
    """
    Sort ``[(start, end)]`` inclusive-inclusive ranges and merge any pair
    that overlaps or is exactly adjacent (``end + 1 == next start``).
    Returns a new list of merged ``(start, end)`` tuples in ascending order.
    """
    if not ranges:
        return []
    ordered = sorted(ranges)
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end + 1:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def memoize(fn):
    """Resettable memoization decorator for instance methods.

    Cached values are stored on each instance instead of in the decorator
    closure. This keeps unrelated instances isolated and prevents the cache
    from extending an instance's lifetime.
    """
    cache_attribute = "_memoize_cache_%s" % fn.__name__

    @wraps(fn)
    def wrapped_fn(self, *args, **kwargs):
        cache = self.__dict__.setdefault(cache_attribute, {})
        cache_key = _memoize_cache_key(args, kwargs)
        try:
            return cache[cache_key]
        except KeyError:
            value = fn(self, *args, **kwargs)
            cache[cache_key] = value
            return value

    def clear_cache(instance):
        instance.__dict__.pop(cache_attribute, None)

    wrapped_fn.clear_cache = clear_cache
    wrapped_fn.make_cache_key = _memoize_cache_key
    return wrapped_fn
