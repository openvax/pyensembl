# TODO: save the config in YMAL file, or TOML file?

MIN_ENSEMBL_RELEASE = 54
MAX_ENSEMBL_RELEASE = 110
MIN_ENSEMBLGENOME_RELEASE = 50
MAX_ENSEMBLGENOME_RELEASE = 57


SPECIES_DATA = [
    {
        "latin_name": "homo_sapiens",
        "synonyms": ["human"],
        "reference_assemblies": {
            "NCBI36": (54, 54),
            "GRCh37": (55, 75),
            "GRCh38": (76, MAX_ENSEMBL_RELEASE),
        },
    },
    {
        "latin_name": "mus_musculus",
        "synonyms": ["mouse", "house mouse"],
        "reference_assemblies": {
            "NCBIM37": (54, 67),
            "GRCm38": (68, 102),
            "GRCm39": (103, MAX_ENSEMBL_RELEASE),
        },
    },
    {
        "latin_name": "canis_familiaris",
        "synonyms": ["dog"],
        "reference_assemblies": {"CanFam3.1": (75, MAX_ENSEMBL_RELEASE)},
    },
    {
        "latin_name": "felis_catus",
        "synonyms": ["cat"],
        "reference_assemblies": {
            "Felis_catus_6.2": (75, 90),
            "Felis_catus_8.0": (91, 92),
            "Felis_catus_9.0": (93, MAX_ENSEMBL_RELEASE),
        },
    },
    {
        "latin_name": "gallus_gallus",
        "synonyms": ["chicken"],
        "reference_assemblies": {
            "Galgal4": (75, 85),
            "Gallus_gallus-5.0": (86, MAX_ENSEMBL_RELEASE),
        },
    },
    {
        "latin_name": "rattus_norvegicus",
        "synonyms": ["rat", "brown_rat", "lab_rat"],
        "reference_assemblies": {
            "Rnor_5.0": (75, 79),
            "Rnor_6.0": (80, 104),
            "mRatBN7.2": (105, MAX_ENSEMBL_RELEASE),
        },
    },
    {
        "latin_name": "macaca_fascicularis",
        "synonyms": ["macaque", "Crab-eating_macaque"],
        "reference_assemblies": {
            "Macaca_fascicularis_6.0": (103, MAX_ENSEMBL_RELEASE)
        },
    },
    {
        "latin_name": "chlorocebus_sabaeus",
        "synonyms": ["green_monkey", "african_green_monkey"],
        "reference_assemblies": {"ChlSab1.1": (86, MAX_ENSEMBL_RELEASE)},
    },
    {
        "latin_name": "macaca_mulatta",
        "synonyms": ["rhesus"],
        "reference_assemblies": {"Mmul_10": (75, MAX_ENSEMBL_RELEASE)},
    },
    {
        "latin_name": "oryctolagus_cuniculus",
        "synonyms": ["rabbit"],
        "reference_assemblies": {"OryCun2.0": (75, MAX_ENSEMBL_RELEASE)},
    },
    {
        "latin_name": "meriones_unguiculatus",
        "synonyms": ["gerbil"],
        "reference_assemblies": {"MunDraft-v1.0": (75, MAX_ENSEMBL_RELEASE)},
    },
    {
        "latin_name": "mesocricetus_auratus",
        "synonyms": ["syrian_hamster"],
        "reference_assemblies": {"MesAur1.0": (75, MAX_ENSEMBL_RELEASE)},
    },
    {
        "latin_name": "cricetulus_griseus_chok1gshd",
        "synonyms": ["chinese_hamster"],
        "reference_assemblies": {"CHOK1GS_HDv1": (75, MAX_ENSEMBL_RELEASE)},
    },
    {
        "latin_name": "heterocephalus_glaber_female",
        "synonyms": ["naked_mole_rat"],
        "reference_assemblies": {
            "HetGla_female_1.0": (75, MAX_ENSEMBL_RELEASE)
        },
    },
    {
        "latin_name": "cavia_porcellus",
        "synonyms": ["guinea_pig"],
        "reference_assemblies": {"Cavpor3.0": (75, MAX_ENSEMBL_RELEASE)},
    },
    {
        "latin_name": "sus_scrofa",
        "synonyms": ["pig"],
        "reference_assemblies": {"Sscrofa11.1": (75, MAX_ENSEMBL_RELEASE)},
    },
    {
        "latin_name": "danio_rerio",
        "synonyms": ["zebrafish"],
        "reference_assemblies": {
            "Zv8": (54, 59),
            "Zv9": (60, 79),
            "GRCz10": (80, 91),
            "GRCz11": (92, MAX_ENSEMBL_RELEASE),
        },
    },
    {
        "latin_name": "drosophila_melanogaster",
        "synonyms": ["drosophila", "fruit fly", "fly"],
        "reference_assemblies": {
            "BDGP5": (75, 78),
            "BDGP6": (79, 95),
            "BDGP6.22": (96, 98),
            "BDGP6.28": (99, 102),
            "BDGP6.32": (103, MAX_ENSEMBL_RELEASE),
        },
    },
    {
        "latin_name": "caenorhabditis_elegans",
        "synonyms": ["nematode", "C_elegans"],
        "reference_assemblies": {
            "WS200": (55, 57),
            "WS210": (58, 60),
            "WS220": (61, 66),
            "WBcel235": (67, MAX_ENSEMBL_RELEASE),
        },
    },
    {
        "latin_name": "saccharomyces_cerevisiae",
        "synonyms": ["yeast", "budding_yeast"],
        "reference_assemblies": {"R64-1-1": (75, MAX_ENSEMBL_RELEASE)},
    },
    {
        "latin_name": "oryza_sativa",
        "synonyms": ["rice", "japanese_rice"],
        "reference_assemblies": {
            "IRGSP-1.0": (55, MAX_ENSEMBLGENOME_RELEASE),
        },
        "database": "plants",
    },
    {
        "latin_name": "arabidopsis_thaliana",
        "synonyms": ["cress", "thale_cress", "hehe"],
        "reference_assemblies": {
            "TAIR10": (55, MAX_ENSEMBLGENOME_RELEASE),
        },
        "database": "plants",
    },
]
