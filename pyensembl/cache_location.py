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

"""Where PyEnsembl keeps its data.

This module only computes the cache root. It is separate from
``download_cache``, which imports hashing, compression and download code, so
that tab completion can find installed releases cheaply.
``pyensembl.download_cache`` re-exports every name defined here.
"""

from os import environ
from os.path import join

CACHE_BASE_SUBDIR = "pyensembl"
CACHE_DIR_ENV_KEY = "PYENSEMBL_CACHE_DIR"


def cache_root():
    """Return the directory holding all of PyEnsembl's data; nothing is created.

    ``$PYENSEMBL_CACHE_DIR/pyensembl`` when that variable is set, otherwise
    the platform cache directory from datacache, e.g. ``~/.cache/pyensembl``,
    ``~/Library/Caches/pyensembl``, or
    ``%LOCALAPPDATA%\\pyensembl\\pyensembl\\Cache``. Every genome's directory
    is ``<root>/<reference>/<annotation><version>``.

    Returns
    -------
    str
        The cache root path.
    """
    environment_root = environ.get(CACHE_DIR_ENV_KEY)
    if environment_root:
        return join(environment_root, CACHE_BASE_SUBDIR)
    import datacache

    return datacache.get_cache_root(CACHE_BASE_SUBDIR)
