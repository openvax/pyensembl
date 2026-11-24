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

from .ensembl_release import EnsemblRelease
from .species import Species, find_species_by_name


def normalize_reference_name(name):
    """
    Search the dictionary of species-specific references (and catalogued
    pinned assemblies) to find a reference name that matches aside from
    capitalization or alias spelling.

    If no matching reference is found, raise an exception.
    """
    lower_name = name.strip().lower()
    for reference in Species.all_reference_names():
        if reference.lower() == lower_name:
            return reference
    from .pinned_genomes import lookup_pinned_alias, resolve_pinned_assembly

    alias = lookup_pinned_alias(name)
    if alias is not None:
        # Raises a detailed ValueError for unsupported versions like T2T-CHM13v1.0.
        return resolve_pinned_assembly(name).reference_name
    raise ValueError("Reference genome '%s' not found" % name)


def find_species_by_reference(reference_name):
    normalized = normalize_reference_name(reference_name)
    species = Species.get_by_reference_name(normalized)
    if species is not None:
        return species
    from .pinned_genomes import resolve_pinned_assembly

    return find_species_by_name(resolve_pinned_assembly(normalized).species)


def which_reference(species_name, ensembl_release):
    return find_species_by_name(species_name).which_reference(ensembl_release)


def genome_for_reference_name(reference_name, allow_older_downloaded_release=True):
    """
    Given a genome reference name, such as "GRCh38" or "T2T-CHM13v2.0"
    (including aliases such as "hs1"), returns the corresponding Ensembl
    release or pinned genome object.

    If `allow_older_downloaded_release` is True, return the newest release
    that is installed (see `Genome.installed`), else the newest whose files
    are downloaded. Choosing only reads the cache; nothing is downloaded.

    Otherwise, or when no release is available locally, return the newest
    release of Ensembl for the reference.
    """
    reference_name = normalize_reference_name(reference_name)
    species = Species.get_by_reference_name(reference_name)
    if species is None:
        from .pinned_genomes import PinnedGenome

        return PinnedGenome.cached(reference_name)
    (min_ensembl_release, max_ensembl_release) = species.reference_assemblies[
        reference_name
    ]
    candidates = [
        EnsemblRelease.cached(release=release, species=species)
        for release in reversed(range(min_ensembl_release, max_ensembl_release + 1))
    ]
    if allow_older_downloaded_release:
        from .release_catalogue import cached_dated_releases
        dated = [g for g in cached_dated_releases()
                 if g.species == species and (g.assembly == reference_name
                                             or g.reference_name == reference_name)]
        candidates += sorted(dated, key=lambda g: g.release, reverse=True)
        for ready in (EnsemblRelease.installed, EnsemblRelease.required_local_files_exist):
            for candidate in candidates:
                if ready(candidate):
                    return candidate
    return candidates[0]


# Convenience genomes and the reference each one selects. ``__getattr__``
# chooses each from the cache the first time it is used, not at import.
CONVENIENCE_GENOME_REFERENCES = {
    "ensembl_grch36": "ncbi36",
    "ensembl_grch37": "grch37",
    "ensembl_grch38": "grch38",
}


def __getattr__(name):
    """Choose a convenience genome the first time it is used (PEP 562).

    ``ensembl_grch36``, ``ensembl_grch37`` and ``ensembl_grch38`` are
    ``genome_for_reference_name`` of NCBI36, GRCh37 and GRCh38, so choosing
    one reads the genome cache. Choosing on first use keeps
    ``import pyensembl`` from touching the cache. The result is stored in the
    module, so later uses return the same object even if the cache changes.

    Parameters
    ----------
    name : str
        ``"ensembl_grch36"``, ``"ensembl_grch37"`` or ``"ensembl_grch38"``.

    Returns
    -------
    EnsemblRelease
        The selected release.

    Raises
    ------
    AttributeError
        For any other name.
    """
    if name not in CONVENIENCE_GENOME_REFERENCES:
        raise AttributeError("module %r has no attribute %r" % (__name__, name))
    genome = genome_for_reference_name(CONVENIENCE_GENOME_REFERENCES[name])
    globals()[name] = genome
    return genome


def __dir__():
    """List the module namespace, including convenience genomes not chosen yet.

    Returns
    -------
    list of str
        Sorted names.
    """
    return sorted(set(globals()) | set(CONVENIENCE_GENOME_REFERENCES))
