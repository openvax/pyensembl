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

from .ensembl_release_versions import MAX_ENSEMBL_RELEASE


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
        for (genome_name, (start, end)) in self.reference_assemblies.items():
            for i in range(start, end + 1):
                assert i not in self._release_to_genome, \
                    "Ensembl release %d already has an associated genome"
                self._release_to_genome[i] = genome_name

    def which_reference(self, ensembl_release):
        if ensembl_release not in self._release_to_genome:
            raise ValueError("No genome for %s in Ensembl release %d" % (
                self.latin_name, ensembl_release))
        return self._release_to_genome[ensembl_release]

    def __str__(self):
        return (
            "Species(latin_name='%s', synonyms=%s, reference_assemblies=%s)" % (
                self.latin_name, self.synonyms, self.reference_assemblies))

    def __repr__(self):
        return str(self)

_latin_names_to_species = {}
_common_names_to_species = {}
_reference_names_to_species = {}

def add_species(latin_name, synonyms, reference_assemblies):
    """
    Create a Species object from the given arguments and enter into
    all the dicts used to look the species up by its fields.
    """
    species = Species(
        latin_name=latin_name,
        synonyms=synonyms,
        reference_assemblies=reference_assemblies)
    _latin_names_to_species[species.latin_name] = species
    for synonym in synonyms:
        if synonym in _common_names_to_species:
            raise ValueError("Can't use synonym '%s' for both %s and %s" % (
                synonym,
                species,
                _common_names_to_species[synonym]))
        _common_names_to_species[synonym] = species
    for reference_name in reference_assemblies:
        if reference_name in _reference_names_to_species:
            raise ValueError("Can't use reference '%s' for both %s and %s" % (
                reference_name,
                species,
                _reference_names_to_species[reference_name]))
        _reference_names_to_species[reference_name] = species
    return species


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

def normalize_reference_name(name):
    """
    Search the dictionary of species-specific references to find a reference
    name that matches aside from capitalization.

    If no matching reference is found, raise an exception.
    """
    lower_name = name.strip().lower()
    for reference in _reference_names_to_species.keys():
        if reference.lower() == lower_name:
            return reference
    raise ValueError("Reference genome '%s' not found" % name)

def find_species_by_name(species_name):
    latin_name = normalize_species_name(species_name)
    if latin_name not in _latin_names_to_species:
        raise ValueError("Species not found: %s" % species_name)
    return _latin_names_to_species[latin_name]

def find_species_by_reference(reference_name):
    return _reference_names_to_species[normalize_reference_name(reference_name)]

def which_reference(species_name, ensembl_release):
    return find_species_by_name(species_name).which_reference(ensembl_release)

def max_ensembl_release(reference_name):
    species = find_species_by_reference(reference_name)
    (_, max_release) = species.reference_assemblies[reference_name]
    return max_release

def check_species_object(species_name_or_object):
    """
    Helper for validating user supplied species names or objects.
    """
    if isinstance(species_name_or_object, Species):
        return species_name_or_object
    elif isinstance(species_name_or_object, str):
        return find_species_by_name(species_name_or_object)
    else:
        raise ValueError("Unexpected type for species: %s : %s" % (
            species_name_or_object, type(species_name_or_object)))

human = add_species(
    latin_name="homo_sapiens",
    synonyms=["human"],
    reference_assemblies={
        "GRCh38": (76, MAX_ENSEMBL_RELEASE),
        "GRCh37": (55, 75),
        "NCBI36": (54, 54),
    })

mouse = add_species(
    latin_name="mus_musculus",
    synonyms=["mouse", "house mouse"],
    # TODO: fix release range
    # The range of releases given for GRCm38  isn't right and we don't have
    # ranges for other assemblies.
    reference_assemblies={
        "GRCm38": (78, MAX_ENSEMBL_RELEASE)
    })

"""
TODO: add all these species

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
