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
    Search the dictionary of species-specific references to find a reference
    name that matches aside from capitalization.

    If no matching reference is found, raise an exception.
    """
    lower_name = name.strip().lower()
    for reference in Species._reference_names_to_species.keys():
        if reference.lower() == lower_name:
            return reference
    raise ValueError("Reference genome '%s' not found" % name)


def find_species_by_reference(reference_name):
    return Species._reference_names_to_species[normalize_reference_name(reference_name)]


def which_reference(species_name, ensembl_release):
    return find_species_by_name(species_name).which_reference(ensembl_release)


def genome_for_reference_name(reference_name, allow_older_downloaded_release=True):
    """
    Given a genome reference name, such as "GRCh38", returns the
    corresponding Ensembl Release object.

    If `allow_older_downloaded_release` is True, return the newest release
    that is installed (see `Genome.installed`), else the newest whose files
    are downloaded. Choosing only reads the cache; nothing is downloaded.

    Otherwise, or when no release is available locally, return the newest
    release of Ensembl for the reference.
    """
    reference_name = normalize_reference_name(reference_name)
    species = find_species_by_reference(reference_name)
    (min_ensembl_release, max_ensembl_release) = species.reference_assemblies[
        reference_name
    ]
    candidates = [
        EnsemblRelease.cached(release=release, species=species)
        for release in reversed(range(min_ensembl_release, max_ensembl_release + 1))
    ]
    if allow_older_downloaded_release:
        for ready in (EnsemblRelease.installed, EnsemblRelease.required_local_files_exist):
            for candidate in candidates:
                if ready(candidate):
                    return candidate
    return candidates[0]


ensembl_grch36 = genome_for_reference_name("ncbi36")
ensembl_grch37 = genome_for_reference_name("grch37")
ensembl_grch38 = genome_for_reference_name("grch38")
