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

from typing import TYPE_CHECKING, List, Literal, Union

from memoized_property import memoized_property

from .locus_with_genome import LocusWithGenome

if TYPE_CHECKING:
    from .exon import Exon
    from .genome import Genome
    from .transcript import Transcript

class Gene(LocusWithGenome):
    def __init__(
        self,
        gene_id: str,
        gene_name: str,
        contig: str,
        start: int,
        end: int,
        strand: Literal["+", "-"],
        biotype: str,
        genome: "Genome",
    ) -> None: ...
    @property
    def id(self) -> str: ...
    @property
    def name(self) -> str: ...
    def __str__(self) -> str: ...
    def __eq__(self, other) -> bool: ...
    def __hash__(self) -> int: ...
    def to_dict(self) -> dict[str, Union[str, int]]: ...
    @memoized_property
    def transcripts(self) -> List[Transcript]: ...
    @memoized_property
    def exons(self) -> list[Exon]: ...
