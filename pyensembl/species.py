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

from __future__ import print_function, absolute_import, division

from .release_info import MIN_ENSEMBL_RELEASE, MAX_ENSEMBL_RELEASE


class Species(object):
    """
    Container for combined information about a species name, its synonyn names
    and which reference to use for this species in each Ensembl release.
    """
    def __init__(self, latin_name, synonyms=[], reference_assemblies={}):
        """
        Parameters
        ----------
        latin_name : str

        synonyms : list of strings

        reference_assemblies : dict
            Mapping of names of reference genomes onto inclusive ranges of
            Ensembl releases Example: {"GRCh37": (54, 75)}
        """
        self.latin_name = latin_name.lower().replace(" ", "_")
        self.synonyms = synonyms
        self.reference_assemblies = reference_assemblies
        self._release_to_genome = {}
        for (genome_name, (start, end)) in self.reference_assemblies:
            for i in range(start, end + 1):
                self._release_to_genome[i] = genome_name

    def which_reference(self, ensembl_release):
        if ensembl_release not in self._release_to_genome:
            raise ValueError("No genome for %s in Ensembl release %d" % (
                self.latin_name, ensembl_release))
        return self._release_to_genome[ensembl_release]


_latin_names_to_species = {}
_common_names_to_species = {}

def add_species(latin_name, synonyms, reference_assemblies):
    species = Species(
        latin_name=latin_name,
        synonyms=synonyms,
        reference_assemblies=reference_assemblies)
    _latin_names_to_species[species.latin_name] = species
    for synonym in synonyms:
        _common_names_to_species[synonym] = species
    return species

grch38 = EnsemblGenomeInfo(
    name="GRCh38",
    first_release=76,
    last_release=MAX_ENSEMBL_RELEASE)

human = add_species(
    latin_name="homo_sapiens",
    synonyms=["human"],
    reference_assemblies=[grch38, grch37, human_ncbi36])

mouse = add_species(
    latin_name="mus_musculus",
    synonyms=["mouse", "house mouse"],
    reference_assemblies=[grcm38])


def normalize_species_name(name):
    """
    If species name was "Homo sapiens" then replace spaces with underscores
    and return "homo_sapiens". Also replace common names like "human" with
    "homo_sapiens".
    """
    lower_name = name.lower().strip()

    # if given a common name such as "human", look up its latin equivalent
    if lower_name in _common_names_to_species:
        return _common_names_to_species[lower_name].latin_name

    return lower_name.replace(" ", "_")

def find_species_by_name(species_name):
    latin_name = normalize_species_name(species_name)
    if latin_name not in _latin_names_to_species:
        raise ValueError("Species not found: %s" % species_name)
    return _latin_names_to_species[latin_name]

def find_species_by_reference(reference_name):
    species = Species._reference_to_species.get(reference_name)
    if not species:
        raise ValueError("Reference genome '%s' not found" % reference_name)
    return species

def species_reference(species_name, ensembl_release):
    return find_species_by_name(species_name).which_reference(ensembl_release)

def max_ensembl_release(reference_name):
    pass

"""

_latin_to_common_names = {
    "ailuropoda_melanoleuca": ["panda"],
    "anas_platyrhynchos": ["mallard", "duck"],
    "anolis_carolinensis": ["carolina anole"],
    "astyanax_mexicanus": ["mexican tetra"],
    "bos_taurus": ["cattle", "cow"],
    "caenorhabditis_elegans": ["c. elegans", "nematode"],
    "callithrix_jacchus": ["marmoset"],
    "canis_familiaris": ["dog"],
    "cavia_porcellus": ["guinea pig"],
    "chlorocebus_sabaeus": ["green monkey", "sabaeus monkey"],
    "choloepus_hoffmanni": ["hoffman's two-toed sloth"],
    "ciona_intestinalis": ["vase tunicate"],
    "ciona_savignyi": [
        "pacific transparent sea squirt",
        "solitary sea squirt"],
    "danio_rerio": ["zebrafish"],
    "dasypus_novemcinctus": ["nine-banded armadillo"],
    "dipodomys_ordii": ["ord's kangaroo rat"],
    "drosophila_melanogaster": [
        "fruit fly",
        "vinegar fly"],
    "echinops_telfairi": ["lesser hedgehog tenrec"],
    "equus_caballus": ["horse"],
    "erinaceus_europaeus": [
        "european hedgehog",
        "hedgehog"],
    "felis_catus": ["cat", "domestic cat"],
    "ficedula_albicollis": ["collared flycatcher"],
    "gadus_morhua": ["atlnatic cod"],
    "gallus_gallus": ["chicken", "domesticated fowl"],
    "gasterosteus_aculeatus": ["three-spined stickleback"],
    "gorilla_gorilla": ["gorilla"],
    human_latin_name: ["human"],
    "ictidomys_tridecemlineatus": [
        "thirteen-lined ground squirrel",
        "striped gopher",
        "leopard ground squirrel"
        "squinney"],
    "latimeria_chalumnae": [
        "west indian ocean coelacanth",
        "african coelacanth"],
    "lepisosteus_oculatus": ["spotted gar"],
    "loxodonta_africana": ["african bush elephant"],
    "macaca_mulatta": ["rhesus macaque"],
    "macropus_eugenii": [
        "tammar wallaby",
        "dama wallaby",
        "darma wallaby"],
    "meleagris_gallopavo": ["wild turkey"],
    "microcebus_murinus": ["gray mouse lemur"],
    "monodelphis_domestica": ["gray short-tailed opossum"],
    mouse_latin_name: ["house mouse", "mouse"],
    "mustela_putorius_furo": ["domestic ferret"],
    "myotis_lucifugus": ["little brown bat"],
    "nomascus_leucogenys": ["northern white-cheeked gibbon"],
    "ochotona_princeps": ["american pika"],
    "oreochromis_niloticus": ["nile tilapia"],
    "ornithorhynchus_anatinus": ["platypus"],
    "oryctolagus_cuniculus": ["european rabbit", "rabbit"],
    "oryzias_latipes": ["japanese rice fish", "medaka"],
    "otolemur_garnettii": [
        "northern greater galago",
        "small-eared greater galago"],
    "ovis_aries": ["sheep"],
    "pan_troglodytes": ["chimpanzee"],
    "papio_anubis": ["olive baboon"],
    "pelodiscus_sinensis": ["chinese softshell turtle"],
    "petromyzon_marinus": ["sea lamprey"],
    "poecilia_formosa": ["amazon molly"],
    "pongo_abelii": ["sumatran orangutan"],
    "procavia_capensis": ["rock hyrax"],
    "pteropus_vampyrus": [
        "large flying fox",
        "large fruit bat",
        "kalang",
        "kalong"],
    "rattus_norvegicus": ["brown rat", "rat", "norwegian rat"],
    "saccharomyces_cerevisiae": ["yeast"],
    "sarcophilus_harrisii": ["tasmanian devil"],
    "sorex_araneus": ["shrew"],
    "sus_scrofa": ["wild boar"],
    "taeniopygia_guttata": ["zebra finch"],
    "takifugu_rubripes": ["japanese puffer", "tiger puffer", "torafugu"],
    "tarsius_syrichta": ["philippine tarsier", "mawmag", "mamag"],
    "tetraodon_nigroviridis": ["green spotted puffer"],
    "tupaia_belangeri": ["northern treeshrew"],
    "tursiops_truncatus": ["bottlenose dolphin"],
    "vicugna_pacos": ["whale"],
    "xenopus_tropicalis": ["western clawed frog"],
    "xiphophorus_maculatus": ["southern platyfish"],
}
"""
