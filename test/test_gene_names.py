"""
Test all methods which return collections of gene names that aren't converting
from some other type of name or ID.
"""
from pyensembl import EnsemblRelease
from test_common import test_ensembl_releases

# make sure that familia
KNOWN_GENE_NAMES = [
    "TP53",
    "ERBB2",
    "SMAD4",
    "CTAG1A",
]

@test_ensembl_releases()
def test_all_gene_names(ensembl):
    """
    test_all_gene_names : Make sure some known gene names such as
    SMAD4, HSP90AA1, TP53, ERBB2
    """
    gene_names = ensembl.gene_names()
    print(type(gene_names))
    for gene_name in KNOWN_GENE_NAMES:
        print(str(gene_name) in set(gene_names))
        assert gene_name in gene_names, \
            "Missing gene name %s from %s" % (gene_name, ensembl)

#
def test_gene_names_at_locus_ensembl77_hla_a():
    # chr6:29,945,884  is a position for HLA-A
    # based on:
    # http://useast.ensembl.org/Homo_sapiens/Gene/
    # Summary?db=core;g=ENSG00000206503;r=6:29941260-29945884
    names = (EnsemblRelease(77, auto_download=True).
             gene_names_at_locus(6, 29945884))
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
