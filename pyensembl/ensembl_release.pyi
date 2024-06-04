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

from typing import Union, TYPE_CHECKING
from typing_extensions import deprecated

from .genome import Genome
from .ensembl_versions import MAX_ENSEMBL_RELEASE
from .species import human
from .ensembl_url_templates import ENSEMBL_FTP_SERVER

if TYPE_CHECKING:
    from .species import Species

class EnsemblRelease(Genome):
    @classmethod
    def normalize_init_values(
        cls, release: Union[int, str], species: Union[Species, str], server: str
    ): ...
    @classmethod
    def cached(
        cls,
        release: int = MAX_ENSEMBL_RELEASE,
        species: Union[str, Species] = human,
        server: str = ENSEMBL_FTP_SERVER,
    ): ...
    def __init__(
        self,
        release: int = MAX_ENSEMBL_RELEASE,
        species: Union[str, Species] = human,
        server: str = ENSEMBL_FTP_SERVER,
    ): ...
    def install_string(self) -> str: ...
    def __str__(self) -> str: ...
    def __eq__(self, other) -> bool: ...
    def __hash__(self) -> int: ...
    def to_dict(self) -> dict: ...
    @classmethod
    def from_dict(cls, state_dict: dict) -> "EnsemblRelease": ...

@deprecated("Use pyensembl.ensembl_release.EnsemblRelease.cached instead.")
def cached_release(release, species="human") -> EnsemblRelease: ...
