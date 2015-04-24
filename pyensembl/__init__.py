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

from .compute_cache import clear_cached_objects
from .ensembl_release import EnsemblRelease
from .gene import Gene
from .gtf import GTF
from .locus import Locus
from .release_info import check_release_number
from .search import find_nearest_locus
from .sequence_data import SequenceData
from .transcript import Transcript

_cache = {}

def cached_release(version, auto_download=True):
    """Cached construction of EnsemblRelease objects. It's desirable to reuse
    the same EnsemblRelease object since each one will store a lot of cached
    annotation data in-memory.
    """
    version = check_release_number(version)
    key = version, auto_download
    if key not in _cache:
        ensembl = EnsemblRelease(version, auto_download=auto_download)
        _cache[key] = ensembl
    return _cache[key]

ensembl_grch36 = ensembl54 = cached_release(54)  # last release for GRCh36/hg18
ensembl_grch37 = ensembl75 = cached_release(75)  # last release for GRCh37/hg19

ensembl77 = cached_release(77)
ensembl78 = cached_release(78)
ensembl79 = cached_release(79)
ensembl_grch38 = ensembl79  # most recent for GRCh38
