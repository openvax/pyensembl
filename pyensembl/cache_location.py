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
SHARED_CACHE_DIR_ENV_KEY = "OPENVAX_DATA_CACHE"

# Directories PyEnsembl itself creates under the cache root which hold no
# genome: shared reference DNA (``genome_fasta_cache``), gene-name alias
# snapshots (``gene_name_aliases``) and dated-release listings
# (``dated_releases``). A scan for installed genomes must skip them.
NON_GENOME_CACHE_DIRECTORIES = frozenset({"dna_cache", "aliases", "dated_releases"})


def resolve_cache_location():
    """Resolve PyEnsembl's cache directory and why it was chosen; nothing is created.

    Precedence:
    1. ``$PYENSEMBL_CACHE_DIR/pyensembl`` when ``PYENSEMBL_CACHE_DIR`` is set.
    2. ``$OPENVAX_DATA_CACHE/pyensembl`` when ``OPENVAX_DATA_CACHE`` is set.
    3. The platform cache directory for ``pyensembl`` from datacache.

    Returns
    -------
    datacache.CacheLocation
        The resolved ``path``, ``origin``, and ``detail``.
    """
    import datacache

    return datacache.resolve_cache_dir(
        CACHE_BASE_SUBDIR,
        CACHE_DIR_ENV_KEY,
        shared_envkey=SHARED_CACHE_DIR_ENV_KEY,
        shared_subdir=None,
        append_package_to_env=True,
    )


def cache_root():
    """Return the directory holding all of PyEnsembl's data; nothing is created.

    ``$PYENSEMBL_CACHE_DIR/pyensembl`` when that variable is set, else
    ``$OPENVAX_DATA_CACHE/pyensembl`` when ``OPENVAX_DATA_CACHE`` is set,
    otherwise the platform cache directory from datacache, e.g.
    ``~/.cache/pyensembl``, ``~/Library/Caches/pyensembl``, or
    ``%LOCALAPPDATA%\\pyensembl\\pyensembl\\Cache``. Every genome's directory
    is ``<root>/<reference>/<annotation><version>``.

    Returns
    -------
    str
        The cache root path.
    """
    environment_root = environ.get(CACHE_DIR_ENV_KEY, "").strip()
    if environment_root:
        from os.path import expanduser

        return join(expanduser(environment_root), CACHE_BASE_SUBDIR)
    shared_root = environ.get(SHARED_CACHE_DIR_ENV_KEY, "").strip()
    if shared_root:
        from os.path import expanduser

        return join(expanduser(shared_root), CACHE_BASE_SUBDIR)
    return resolve_cache_location().path
