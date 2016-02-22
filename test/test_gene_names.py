"""
Test all methods which return collections of gene names that aren't converting
from some other type of name or ID.
"""
from __future__ import absolute_import, print_function
from pyensembl import ensembl_grch38

from .common import test_ensembl_releases

KNOWN_GENE_NAMES = [
    "TP53",
    "ERBB2",
    "SMAD4",
    "CTAG1A",
    "HLA-A",
]

@test_ensembl_releases()
def test_all_gene_names(ensembl):
    """
    test_all_gene_names : Make sure some known gene names such as
    SMAD4, TP53, ERBB2, &c
    """
    gene_names = ensembl.gene_names()
    print(type(gene_names))
    for gene_name in KNOWN_GENE_NAMES:
        assert gene_name in gene_names, \
            "Missing gene name %s from %s" % (gene_name, ensembl)

def test_gene_names_at_locus_grch38_hla_a():
    # chr6:29,945,884  is a position for HLA-A
    # based on:
    # http://useast.ensembl.org/Homo_sapiens/Gene/
    # Summary?db=core;g=ENSG00000206503;r=6:29941260-29945884
    names = ensembl_grch38.gene_names_at_locus(6, 29945884)
    assert names == ["HLA-A"], "Expected gene name HLA-A, got: %s" % (names,)

@test_ensembl_releases()
def test_gene_names_on_contig(ensembl):
    gene_names_chr17 = ensembl.gene_names(17)
    assert "TP53" in gene_names_chr17, \
        "No TP53 in gene names on chr17 of %s, gene names: %s ... (%d)" % (
            ensembl, list(gene_names_chr17[:4]), len(gene_names_chr17))

    gene_names_chr18 = ensembl.gene_names(18)
    assert "SMAD4" in gene_names_chr18, \
        "No SMAD4 in gene names on chr18 of %s, gene names: %s ... (%d)" % (
            ensembl, list(gene_names_chr18[:4]), len(gene_names_chr18))


def test_gene_name_of_HLA_gene_id():
    gene_ids = ensembl_grch38.gene_ids_of_gene_name("HLA-A")
    gene_names = [
        ensembl_grch38.gene_name_of_gene_id(gene_id)
        for gene_id in gene_ids
    ]
    unique_gene_names = list(set(gene_names))
    assert len(unique_gene_names) == 1, (len(unique_gene_names), unique_gene_names)
    gene_name = unique_gene_names[0]
    assert gene_name == "HLA-A", gene_name
