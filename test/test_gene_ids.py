"""
Test all methods which return collections of gene IDs that aren't converting
from some other tyoe of name or ID.

TODO: Implement tests for EnsemblRelease.gene_ids
"""
from __future__ import absolute_import
from .test_common import test_ensembl_releases
from pyensembl import ensembl_grch38 as ensembl

def test_gene_ids_ensembl77_hla_a():
    # chr6:29,945,884  is a position for HLA-A
    # Gene ID = ENSG00000206503
    # based on:
    # http://useast.ensembl.org/Homo_sapiens/Gene/
    # Summary?db=core;g=ENSG00000206503;r=6:29941260-29945884
    ids = ensembl.gene_ids_at_locus(6, 29945884)
    expected = "ENSG00000206503"
    assert ids == ["ENSG00000206503"], \
        "Expected HLA-A, gene ID = %s, got: %s" % (expected, ids)

@test_ensembl_releases()
def test_gene_ids_on_contig(ensembl):
    gene_ids_chr17 = ensembl.gene_ids(contig=17)
    # gene ID of TP53
    tp53 = "ENSG00000141510"
    assert tp53 in gene_ids_chr17, \
        "Missing %s from %s on chr17, example IDs: %s (total = %d)" % (
            tp53, ensembl, gene_ids_chr17[:5], len(gene_ids_chr17))

    # gene ID of SMAD4
    gene_ids_chr18 = ensembl.gene_ids(contig=18)
    smad4 = "ENSG00000141646"
    assert smad4 in gene_ids_chr18, \
        "Missing %s from %s on chr18, example result: %s (total = %d)" % (
            smad4, ensembl, gene_ids_chr18[:5], len(gene_ids_chr18))

