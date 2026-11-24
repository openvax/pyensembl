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

"""Implementation of the ``pyensembl install`` CLI action."""

import logging

from .formatting import genome_description

logger = logging.getLogger(__name__)


class InstallError(Exception):
    """An install failed because a source file is missing."""


def progress_available():
    """Return whether CLI progress bars should be displayed.

    Delegates to ``datacache.should_show_progress("auto")``, which enables
    progress bars when stderr is an interactive TTY and neither
    ``OPENVAX_NO_PROGRESS`` nor ``DATACACHE_NO_PROGRESS`` disables them.

    Returns
    -------
    bool
        ``True`` if progress bars should be shown, ``False`` otherwise.
    """
    import datacache

    return datacache.should_show_progress("auto")


def install_genome(genome, only_genome_fasta=False, overwrite=False, show_progress=False):
    """Download and index a single genome (or only its reference DNA).

    Parameters
    ----------
    genome : Genome
        Genome instance to install.
    only_genome_fasta : bool, optional
        When ``True``, download and index only reference DNA.
    overwrite : bool, optional
        When ``True``, redo downloads and indexing even if files exist.
    show_progress : bool, optional
        When ``True``, show progress bars during download and indexing.
    """
    description = genome_description(genome)
    if only_genome_fasta:
        installed = genome.genome_fasta_status() == "indexed"
    else:
        installed = genome.installed()
    if installed and not overwrite:
        # The steps below then only confirm that everything is in place.
        logger.info("%s is already installed", description)
    else:
        logger.info("Installing %s", description)
    if only_genome_fasta:
        genome.download_genome_fasta(overwrite=overwrite, show_progress=show_progress)
        genome.index_genome_fasta(overwrite=overwrite, show_progress=show_progress)
    else:
        genome.download(overwrite=overwrite, show_progress=show_progress)
        genome.index(overwrite=overwrite, show_progress=show_progress)


def run_install(genome, only_genome_fasta=False, overwrite=False):
    """Download and index one genome, reporting missing source files as ``InstallError``.

    Parameters
    ----------
    genome : Genome
        Genome to download and index.
    only_genome_fasta : bool, optional
        Install only reference DNA.
    overwrite : bool, optional
        Redo downloads and indexing even when files exist.

    Raises
    ------
    UnpublishedDateError
        If a dated release is not published.
    InstallError
        For missing local or remote source files; ``str(error)`` is the
        user-facing error message.
    """
    from ..download_cache import MissingGenomeDataFile, MissingLocalFile
    from ..genome_fasta import MissingGenomeFastaError

    try:
        install_genome(
            genome,
            only_genome_fasta,
            overwrite,
            show_progress=progress_available(),
        )
    except (MissingLocalFile, MissingGenomeDataFile, MissingGenomeFastaError) as error:
        # Missing remote cache entries are converted to the narrow
        # MissingGenomeDataFile type by DownloadCache. Parser,
        # integrity and database failures keep their tracebacks.
        message = (
            "file not found: %s" % error.path
            if isinstance(error, MissingLocalFile)
            else str(error)
        )
        raise InstallError(message) from error
