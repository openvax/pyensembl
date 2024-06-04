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

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .ensembl_release import EnsemblRelease
    from .species import Species

def normalize_reference_name(name: str) -> str: ...
def find_species_by_reference(reference_name: str) -> Species: ...
def which_reference(species_name: str, ensembl_release: int) -> str: ...
def max_ensembl_release(reference_name: str) -> int: ...
def genome_for_reference_name(
    reference_name: str, allow_older_downloaded_release: bool = True
) -> EnsemblRelease: ...

ensembl_grch36: EnsemblRelease = ...
ensembl_grch37: EnsemblRelease = ...
ensembl_grch38: EnsemblRelease = ...
