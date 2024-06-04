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

from typing import Dict, Generator, List, Tuple, Union

from serializable import Serializable

from .ensembl_versions import MAX_ENSEMBL_RELEASE, MAX_PLANTS_ENSEMBL_RELEASE

class Species(Serializable):
    _latin_names_to_species: dict[str, Species] = {}
    _common_names_to_species: dict[str, Species] = {}
    _reference_names_to_species: dict[str, Species] = {}

    @classmethod
    def register(
        cls,
        latin_name: str,
        synonyms: List[str],
        reference_assemblies: Dict[str, Tuple[int, int]],
        is_plant: bool = False,
    ) -> Species: ...
    @classmethod
    def all_registered_latin_names(cls) -> List[str]: ...
    @classmethod
    def all_species_release_pairs(cls) -> Generator[Tuple[str, int], None, None]: ...
    def __init__(
        self,
        latin_name: str,
        synonyms: list[str] = [],
        reference_assemblies: dict[str, tuple[int, int]] = {},
        is_plant: bool = False,
    ): ...
    def which_reference(self, ensembl_release: int) -> str: ...
    def __str__(self) -> str: ...
    def __eq__(self, other) -> bool: ...
    def to_dict(self) -> Dict[str, str]: ...
    @classmethod
    def from_dict(cls, state_dict: Dict[str, str]) -> Species: ...
    def __hash__(self) -> int: ...

def normalize_species_name(name: str) -> str: ...
def find_species_by_name(species_name: str) -> Species: ...
def check_species_object(species_name_or_object: Union[str, Species]) -> Species: ...

human = Species.register(
    latin_name="homo_sapiens",
    synonyms=["human"],
    reference_assemblies={
        "GRCh38": (76, MAX_ENSEMBL_RELEASE),
        "GRCh37": (55, 75),
        "NCBI36": (54, 54),
    },
)

mouse = Species.register(
    latin_name="mus_musculus",
    synonyms=["mouse", "house mouse"],
    reference_assemblies={
        "NCBIM37": (54, 67),
        "GRCm38": (68, 102),
        "GRCm39": (103, MAX_ENSEMBL_RELEASE),
    },
)

dog = Species.register(
    latin_name="canis_familiaris",
    synonyms=["dog"],
    reference_assemblies={"CanFam3.1": (75, MAX_ENSEMBL_RELEASE)},
)

cat = Species.register(
    latin_name="felis_catus",
    synonyms=["cat"],
    reference_assemblies={
        "Felis_catus_6.2": (75, 90),
        "Felis_catus_8.0": (91, 92),
        "Felis_catus_9.0": (93, MAX_ENSEMBL_RELEASE),
    },
)

chicken = Species.register(
    latin_name="gallus_gallus",
    synonyms=["chicken"],
    reference_assemblies={
        "Galgal4": (75, 85),
        "Gallus_gallus-5.0": (86, MAX_ENSEMBL_RELEASE),
    },
)

# Does the black rat (Rattus Rattus) get used for research too?
brown_rat = Species.register(
    latin_name="rattus_norvegicus",
    synonyms=["brown rat", "lab rat", "rat"],
    reference_assemblies={
        "Rnor_5.0": (75, 79),
        "Rnor_6.0": (80, 104),
        "mRatBN7.2": (105, MAX_ENSEMBL_RELEASE),
    },
)

macaque = Species.register(
    latin_name="macaca_fascicularis",
    synonyms=["macaque", "Crab-eating macaque"],
    reference_assemblies={
        "Macaca_fascicularis_6.0": (103, MAX_ENSEMBL_RELEASE),
    },
)

green_monkey = Species.register(
    latin_name="chlorocebus_sabaeus",
    synonyms=["green_monkey", "african_green_monkey"],
    reference_assemblies={
        "ChlSab1.1": (86, MAX_ENSEMBL_RELEASE),
    },
)

rhesus = Species.register(
    latin_name="macaca_mulatta",
    synonyms=["rhesus"],
    reference_assemblies={"Mmul_10": (75, MAX_ENSEMBL_RELEASE)},
)

rabbit = Species.register(
    latin_name="oryctolagus_cuniculus",
    synonyms=["rabbit"],
    reference_assemblies={"OryCun2.0": (75, MAX_ENSEMBL_RELEASE)},
)

gerbil = Species.register(
    latin_name="meriones_unguiculatus",
    synonyms=["gerbil"],
    reference_assemblies={"MunDraft-v1.0": (75, MAX_ENSEMBL_RELEASE)},
)

syrian_hamster = Species.register(
    latin_name="mesocricetus_auratus",
    synonyms=["syrian_hamster"],
    reference_assemblies={"MesAur1.0": (75, MAX_ENSEMBL_RELEASE)},
)

chinese_hamster = Species.register(
    latin_name="cricetulus_griseus_chok1gshd",
    synonyms=["chinese_hamster"],
    reference_assemblies={"CHOK1GS_HDv1": (75, MAX_ENSEMBL_RELEASE)},
)

naked_mole_rat = Species.register(
    latin_name="heterocephalus_glaber_female",
    synonyms=["naked_mole_rat"],
    reference_assemblies={"HetGla_female_1.0": (75, MAX_ENSEMBL_RELEASE)},
)

guinea_pig = Species.register(
    latin_name="cavia_porcellus",
    synonyms=["guinea_pig"],
    reference_assemblies={"Cavpor3.0": (75, MAX_ENSEMBL_RELEASE)},
)

pig = Species.register(
    latin_name="sus_scrofa",
    synonyms=["pig"],
    reference_assemblies={"Sscrofa11.1": (75, MAX_ENSEMBL_RELEASE)},
)

zebrafish = Species.register(
    latin_name="danio_rerio",
    synonyms=["zebrafish"],
    reference_assemblies={
        "ZFISH7": (47, 53),
        "Zv8": (54, 59),
        "Zv9": (60, 79),
        "GRCz10": (80, 91),
        "GRCz11": (92, MAX_ENSEMBL_RELEASE),
    },
)

fly = Species.register(
    latin_name="drosophila_melanogaster",
    synonyms=["drosophila", "fruit fly", "fly"],
    reference_assemblies={
        "BDGP5": (75, 78),
        "BDGP6": (79, 95),
        "BDGP6.22": (96, 98),
        "BDGP6.28": (99, 102),
        "BDGP6.32": (103, MAX_ENSEMBL_RELEASE),
    },
)

nematode = Species.register(
    latin_name="caenorhabditis_elegans",
    synonyms=["nematode", "C_elegans"],
    reference_assemblies={
        "WS180": (47, 49),
        "WS190": (50, 54),
        "WS200": (55, 57),
        "WS210": (58, 59),
        "WS220": (61, 66),
        "WBcel215": (67, 70),
        "WBcel235": (71, MAX_ENSEMBL_RELEASE),
    },
)

yeast = Species.register(
    latin_name="saccharomyces_cerevisiae",
    synonyms=["yeast", "budding_yeast"],
    reference_assemblies={
        "R64-1-1": (76, MAX_ENSEMBL_RELEASE),
    },
)

arabidopsis_thaliana = Species.register(
    latin_name="arabidopsis_thaliana",
    synonyms=["arabidopsis"],
    reference_assemblies={
        "TAIR10": (40, MAX_PLANTS_ENSEMBL_RELEASE),
    },
    is_plant=True,
)

rice = Species.register(
    latin_name="oryza_sativa",
    synonyms=["rice"],
    reference_assemblies={
        "IRGSP-1.0": (40, MAX_PLANTS_ENSEMBL_RELEASE),
    },
    is_plant=True,
)
