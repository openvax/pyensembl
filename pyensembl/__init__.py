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

from __future__ import print_function, division, absolute_import

from .memory_cache import MemoryCache
from .download_cache import DownloadCache
from .ensembl_release import EnsemblRelease
from .ensembl_release_versions import check_release_number, MAX_ENSEMBL_RELEASE
from .genome import Genome
from .gene import Gene
from .gtf import GTF
from .locus import Locus
from .search import find_nearest_locus
from .sequence_data import SequenceData
from .species import (
    find_species_by_name,
    find_species_by_reference,
    which_reference,
    check_species_object,
    normalize_reference_name,
    normalize_species_name,
)
from .transcript import Transcript

_cache = {}

def cached_release(version, species="human"):
    """Cached construction of EnsemblRelease objects. It's desirable to reuse
    the same EnsemblRelease object since each one will store a lot of cached
    annotation data in-memory.
    """
    version = check_release_number(version)
    species = check_species_object(species)
    key = (version, species)
    if key not in _cache:
        ensembl = EnsemblRelease(version, species=species)
        _cache[key] = ensembl
    return _cache[key]

def genome_for_reference_name(reference_name):
    reference_name = normalize_reference_name(reference_name)
    species = find_species_by_reference(reference_name)
    (_, max_ensembl_release) = species.reference_assemblies[reference_name]
    return cached_release(max_ensembl_release, species=species)

ensembl_grch36 = ensembl54 = cached_release(54)  # last release for GRCh36/hg18
ensembl_grch37 = ensembl75 = cached_release(75)  # last release for GRCh37/hg19

ensembl77 = cached_release(77)
ensembl78 = cached_release(78)
ensembl79 = cached_release(79)
ensembl80 = cached_release(80)
ensembl81 = cached_release(81)
ensembl_grch38 = ensembl81  # most recent for GRCh38
