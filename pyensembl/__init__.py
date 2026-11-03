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

"""Python interface to Ensembl reference genome metadata.

Public names are imported from their submodules on first access (PEP 562),
so ``import pyensembl`` stays fast and does not read the genome cache.
``from pyensembl import EnsemblRelease`` and ``pyensembl.EnsemblRelease``
work as before.
"""

import importlib
import logging
from typing import TYPE_CHECKING

from .version import __version__

if TYPE_CHECKING:
    from .database import Database
    from .dated_releases import (
        UnpublishedDateError,
        available_annotation_dates,
        available_dated_releases,
    )
    from .release_catalogue import available_releases
    from .download_cache import DownloadCache
    from .ensembl_release import EnsemblRelease, cached_release
    from .ensembl_versions import MAX_ENSEMBL_RELEASE
    from .exon import Exon
    from .genome import Genome
    from .genome_fasta import MissingGenomeFastaError
    from .genome_fasta_cache import prune_genome_fastas
    from .gene import Gene
    from .genetic_code import GeneticCodeAssignment
    from .gene_name_aliases import GeneNameAliases
    from .locus import Locus
    from .reference_name import (
        ensembl_grch36,
        ensembl_grch37,
        ensembl_grch38,
        normalize_reference_name,
        find_species_by_reference,
        which_reference,
        genome_for_reference_name,
    )
    from .protein import Protein
    from .search import find_nearest_locus
    from .sequence_data import SequenceData
    from .species import find_species_by_name, check_species_object, normalize_species_name
    from .transcript import Transcript

# Per the Python logging HOWTO ("Configuring Logging for a Library"), attach a
# no-op handler to the package logger so that library usage neither emits log
# output nor triggers "No handlers could be found" warnings, while leaving the
# root logger and any application-configured loggers untouched.
logging.getLogger(__name__).addHandler(logging.NullHandler())

# The submodule that defines each public name. ``__getattr__`` imports it the
# first time the name is used.
EXPORT_MODULES = {
    "Database": "database",
    "UnpublishedDateError": "dated_releases",
    "available_annotation_dates": "dated_releases",
    "available_dated_releases": "dated_releases",
    "available_releases": "release_catalogue",
    "DownloadCache": "download_cache",
    "EnsemblRelease": "ensembl_release",
    "cached_release": "ensembl_release",
    "MAX_ENSEMBL_RELEASE": "ensembl_versions",
    "Exon": "exon",
    "Genome": "genome",
    "MissingGenomeFastaError": "genome_fasta",
    "prune_genome_fastas": "genome_fasta_cache",
    "Gene": "gene",
    "GeneticCodeAssignment": "genetic_code",
    "GeneNameAliases": "gene_name_aliases",
    "Locus": "locus",
    # Chosen from the genome cache when first used; see reference_name.py.
    "ensembl_grch36": "reference_name",
    "ensembl_grch37": "reference_name",
    "ensembl_grch38": "reference_name",
    "normalize_reference_name": "reference_name",
    "find_species_by_reference": "reference_name",
    "which_reference": "reference_name",
    "genome_for_reference_name": "reference_name",
    "Protein": "protein",
    "find_nearest_locus": "search",
    "SequenceData": "sequence_data",
    "find_species_by_name": "species",
    "check_species_object": "species",
    "normalize_species_name": "species",
    "Transcript": "transcript",
}

__all__ = [
    "__version__",
    "available_releases",
    "available_annotation_dates",
    "available_dated_releases",
    "UnpublishedDateError",
    "DownloadCache",
    "Database",
    "EnsemblRelease",
    "cached_release",
    "MAX_ENSEMBL_RELEASE",
    "Gene",
    "GeneticCodeAssignment",
    "GeneNameAliases",
    "Transcript",
    "Protein",
    "Exon",
    "SequenceData",
    "find_nearest_locus",
    "find_species_by_name",
    "find_species_by_reference",
    "genome_for_reference_name",
    "which_reference",
    "check_species_object",
    "normalize_reference_name",
    "normalize_species_name",
    "Genome",
    "MissingGenomeFastaError",
    "prune_genome_fastas",
    "Locus",
    "Exon",
    "ensembl_grch36",
    "ensembl_grch37",
    "ensembl_grch38",
]


def __getattr__(name):
    """Import a public name or a submodule the first time it is used (PEP 562).

    Python calls this only for names not yet in the package namespace. The
    result is stored there, so later lookups do not come back here.

    Parameters
    ----------
    name : str
        A name from ``__all__`` or the name of a submodule, e.g. ``"shell"``.

    Returns
    -------
    object
        The public object, or the submodule ``pyensembl.<name>``.

    Raises
    ------
    AttributeError
        If ``name`` is neither a public name nor a submodule.
    """
    if name in EXPORT_MODULES:
        module = importlib.import_module("." + EXPORT_MODULES[name], __name__)
        value = getattr(module, name)
    elif name.startswith("_"):
        raise AttributeError("module %r has no attribute %r" % (__name__, name))
    else:
        # Before 2.33, ``import pyensembl`` imported most submodules, so code
        # could use e.g. ``pyensembl.species`` without importing it.
        try:
            value = importlib.import_module("." + name, __name__)
        except ModuleNotFoundError as error:
            if error.name != "%s.%s" % (__name__, name):
                raise
            raise AttributeError(
                "module %r has no attribute %r" % (__name__, name)
            ) from None
    globals()[name] = value
    return value


def __dir__():
    """List the package namespace, including public names not imported yet.

    Returns
    -------
    list of str
        Sorted names, as ``dir(pyensembl)`` and interactive completion show them.
    """
    return sorted(set(globals()) | set(__all__))
