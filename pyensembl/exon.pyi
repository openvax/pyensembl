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

from typing import Dict, Literal, Union
from .locus import Locus

class Exon(Locus):
    def __init__(
        self,
        exon_id: str,
        contig: str,
        start: int,
        end: int,
        strand: Literal["+", "-"],
        gene_name: str,
        gene_id: str,
    ) -> None: ...
    @property
    def id(self) -> str: ...
    def __str__(self) -> str: ...
    def __eq__(self, other) -> bool: ...
    def __hash__(self) -> int: ...
    def to_dict(self) -> Dict[str, Union[str, int]]: ...
