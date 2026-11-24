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

"""Implementation of ``delete-index-files``, ``delete-all-files``, and ``prune``."""

import os

from .formatting import directory_size, genome_description


def delete_genome_files(genome, action):
    """Delete index files or all cached files for ``genome`` and print what was removed.

    Parameters
    ----------
    genome : Genome
        Target genome whose files should be deleted.
    action : str
        Either ``"delete-index-files"`` or ``"delete-all-files"``.
    """
    from ..genome_fasta_cache import dna_cache_lock

    if action == "delete-index-files":
        deleted = genome.delete_index_files()
    else:
        directory = genome.download_cache.cache_directory_path
        deleted = []
        # Serialize reference removal with shared-DNA install/prune. Shared
        # objects themselves are removed only by the explicit prune action.
        if os.path.isdir(directory):
            with dna_cache_lock():
                if os.path.isdir(directory):
                    genome.close()
                    size = directory_size(directory)
                    genome.download_cache.delete_cache_directory()
                    deleted.append((directory, size))
    for path, size in deleted:
        print("Deleted %s (%s bytes)" % (path, format(size, ",")))
    if not deleted:
        print("Nothing to delete for %s" % genome_description(genome))


def run_prune(dry_run):
    """Delete (or with ``dry_run`` only list) reference DNA no release uses.

    Parameters
    ----------
    dry_run : bool
        Report candidates without deleting files.

    Returns
    -------
    str
        One line per candidate, or a note that there are none.

    Raises
    ------
    ValueError
        If the shared DNA cache cannot be pruned safely.
    """
    from ..genome_fasta_cache import prune_genome_fastas

    candidates = prune_genome_fastas(dry_run=dry_run)
    if not candidates:
        return "No orphan genome FASTAs"
    verb = "Would delete" if dry_run else "Deleted"
    return "\n".join(
        "%s %s (%s bytes)" % (verb, path, format(size, ",")) for path, size in candidates
    )
