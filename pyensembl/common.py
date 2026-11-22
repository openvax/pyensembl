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
import hashlib
import os
from pathlib import Path
import pickle
import re
from uuid import uuid4


_MD5SUM_LINE = re.compile(r"([a-fA-F0-9]{32})\s+\*?(?:\./)?([^/]+)")


def _md5_hasher():
    return hashlib.md5(usedforsecurity=False)


def file_digest(path, algorithm="sha256", *, opener=open):
    """Return the hex digest of a file's bytes."""
    digest = _md5_hasher if algorithm == "md5" else algorithm
    with opener(path, "rb") as handle:
        return hashlib.file_digest(handle, digest).hexdigest()


def parse_md5sum_manifest(text, filename):
    """Return the unique lowercase MD5 hex digest for ``filename`` in ``text``."""
    matches = []
    for line in text.splitlines():
        match = _MD5SUM_LINE.fullmatch(line.strip())
        if match and match[2] == filename:
            matches.append(match[1].lower())
    if len(matches) != 1:
        raise ValueError("Expected one checksum for %s" % filename)
    return matches[0]


def file_fingerprint(path, *, ctime=True, include_dev=False):
    """Return JSON-compatible stat identity for a file.

    Returns ``[size, mtime_ns, ctime_ns, ino]`` by default, omitting
    ``ctime_ns`` when ``ctime=False`` (for cache-owned files whose metadata
    may change on chmod/chown) and appending ``st_dev`` when
    ``include_dev=True``.
    """
    info = os.stat(path)
    fields = [info.st_size, info.st_mtime_ns]
    if ctime:
        fields.append(info.st_ctime_ns)
    fields.append(info.st_ino)
    if include_dev:
        fields.append(info.st_dev)
    return fields


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


@contextmanager
def _read_progress(raw, description, enabled):
    """Yield a callable that reports how much of an open binary file was read."""
    if not enabled:
        yield lambda: None
        return
    from tqdm.auto import tqdm

    total = os.fstat(raw.fileno()).st_size
    with tqdm(
        total=total, desc=description, unit="B", unit_scale=True, unit_divisor=1024,
        leave=False,
    ) as bar:
        yield lambda: bar.update(raw.tell() - bar.n)


def dump_pickle(obj, filepath):
    # Atomic, so an interrupted write never leaves a truncated index behind.
    with _atomic_output(filepath) as f:
        pickle.dump(obj, file=f)


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
