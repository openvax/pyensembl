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

_latin_to_common_names = {
    "ailuropoda_melanoleuca" : ["panda"],
    "anas_platyrhynchos" : ["mallard", "duck"],
    "anolis_carolinensis" : ["carolina anole"],
    "astyanax_mexicanus" : ["mexican tetra"],
    "bos_taurus" : ["cattle", "cow"],
    "caenorhabditis_elegans" : ["c. elegans", "nematode"],
    "callithrix_jacchus" : ["marmoset"],
    "canis_familiaris" : ["dog"],
    "cavia_porcellus" : ["guinea pig"],
    "chlorocebus_sabaeus" : ["green monkey", "sabaeus monkey"],
    "choloepus_hoffmanni" : ["hoffman's two-toed sloth"],
    "ciona_intestinalis" : ["vase tunicate"],
    "ciona_savignyi" : [
        "pacific transparent sea squirt",
        "solitary sea squirt"],
    "danio_rerio" : ["zebrafish"],
    "dasypus_novemcinctus" : ["nine-banded armadillo"],
    "dipodomys_ordii" : ["ord's kangaroo rat"],
    "drosophila_melanogaster" : [
        "fruit fly",
        "vinegar fly"],
    "echinops_telfairi" : ["lesser hedgehog tenrec"],
    "equus_caballus" : ["horse"],
    "erinaceus_europaeus" : [
        "european hedgehog",
        "hedgehog"],
    "felis_catus" : ["cat", "domestic cat"],
    "ficedula_albicollis" : ["collared flycatcher"],
    "gadus_morhua" : ["atlnatic cod"],
    "gallus_gallus" : ["chicken", "domesticated fowl"],
    "gasterosteus_aculeatus" : ["three-spined stickleback"],
    "gorilla_gorilla" : ["gorilla"],
    "homo_sapiens" : ["human"],
    "ictidomys_tridecemlineatus" : [
        "thirteen-lined ground squirrel",
        "striped gopher",
        "leopard ground squirrel"
        "squinney"],
    "latimeria_chalumnae" : [
        "west indian ocean coelacanth",
        "african coelacanth"],
    "lepisosteus_oculatus" : ["spotted gar"],
    "loxodonta_africana" : ["african bush elephant"],
    "macaca_mulatta" : ["rhesus macaque"],
    "macropus_eugenii" : [
        "tammar wallaby",
        "dama wallaby",
        "darma wallaby"],
    "meleagris_gallopavo" : ["wild turkey"],
    "microcebus_murinus" : ["gray mouse lemur"],
    "monodelphis_domestica" : ["gray short-tailed opossum"],
    "mus_musculus" : ["house mouse", ],
    "mustela_putorius_furo" : ["domestic ferret"],
    "myotis_lucifugus" : ["little brown bat"],
    "nomascus_leucogenys" : ["northern white-cheeked gibbon"],
    "ochotona_princeps" : ["american pika"],
    "oreochromis_niloticus" : ["nile tilapia"],
    "ornithorhynchus_anatinus" : ["platypus"],
    "oryctolagus_cuniculus" : ["european rabbit", "rabbit"],
    "oryzias_latipes" : ["japanese rice fish", "medaka"],
    "otolemur_garnettii" : [
        "northern greater galago",
        "small-eared greater galago"],
    "ovis_aries" : ["sheep"],
    "pan_troglodytes" : ["chimpanzee"],
    "papio_anubis" : ["olive baboon"],
    "pelodiscus_sinensis" : ["chinese softshell turtle"],
    "petromyzon_marinus" : ["sea lamprey"],
    "poecilia_formosa" : ["amazon molly"],
    "pongo_abelii" : ["sumatran orangutan"],
    "procavia_capensis" : ["rock hyrax"],
    "pteropus_vampyrus" : [
        "large flying fox",
        "large fruit bat",
        "kalang",
        "kalong"],
    "rattus_norvegicus" : ["brown rat", "rat", "norwegian rat"],
    "saccharomyces_cerevisiae" : ["yeast"],
    "sarcophilus_harrisii" : ["tasmanian devil"],
    "sorex_araneus" : ["shrew"],
    "sus_scrofa" : ["wild boar"],
    "taeniopygia_guttata" : ["zebra finch"],
    "takifugu_rubripes" : ["japanese puffer", "tiger puffer", "torafugu"],
    "tarsius_syrichta" : ["philippine tarsier", "mawmag", "mamag"],
    "tetraodon_nigroviridis" : ["green spotted puffer"],
    "tupaia_belangeri" : ["northern treeshrew"],
    "tursiops_truncatus" : ["bottlenose dolphin"],
    "vicugna_pacos" : ["whale"],
    "xenopus_tropicalis" : ["western clawed frog"],
    "xiphophorus_maculatus" : ["southern platyfish"],
}

_common_to_latin_name = {}
for (latin_name, common_names) in _latin_to_common_names.items():
    for common_name in common_names:
        assert common_name not in _common_to_latin_name, \
            "Common name %s (for %s) appears twice" % (common_name, latin_name)
        _common_to_latin_name[common_name] = latin_name

def normalize_species_name(name):
    name = name.lower()

    # if name is already in form like "homo_sapiens" then just return it
    if name in _latin_to_common_names:
        return name

    # if species name was "homo sapiens" then replace spaces with underscores
    # and return "homo_sapiens"
    if " " in name:
        no_spaces = name.replace(" ", "_")
        if no_spaces in _latin_to_common_names:
            return no_spaces

    # if given a common name such as "human", look up its latin equivalent
    if name in _common_to_latin_name:
        return _common_to_latin_name[name]

    raise ValueError("Unknown species name: %s" % name)

