"""
Test all methods which return collections of gene IDs that aren't converting
from some other type of name or ID.

TODO: Implement tests for EnsemblRelease.gene_ids
"""
from __future__ import absolute_import

from nose.tools import assert_raises, ok_
from pyensembl import ensembl_grch38, cached_release

from .common import test_ensembl_releases

ensembl77 = cached_release(77, "human")

def test_gene_ids_grch38_hla_a():
    # chr6:29,945,884  is a position for HLA-A
    # Gene ID = ENSG00000206503
    # based on:
    # http://useast.ensembl.org/Homo_sapiens/Gene/
    # Summary?db=core;g=ENSG00000206503;r=6:29941260-29945884
    ids = ensembl_grch38.gene_ids_at_locus(6, 29945884)
    expected = "ENSG00000206503"
    assert ids == ["ENSG00000206503"], \
        "Expected HLA-A, gene ID = %s, got: %s" % (expected, ids)

def test_gene_ids_of_gene_name_hla_grch38():
    hla_a_gene_ids = ensembl_grch38.gene_ids_of_gene_name("HLA-A")
    assert 'ENSG00000206503' in hla_a_gene_ids, hla_a_gene_ids

    hla_b_gene_ids = ensembl_grch38.gene_ids_of_gene_name("HLA-B")
    assert 'ENSG00000234745' in hla_b_gene_ids, hla_b_gene_ids

    hla_c_gene_ids = ensembl_grch38.gene_ids_of_gene_name("HLA-C")
    assert 'ENSG00000204525' in hla_c_gene_ids, hla_c_gene_ids

def test_gene_id_of_protein_id_release77():
    gene_id = ensembl77.gene_id_of_protein_id("ENSP00000485677")
    ok_('ENSG00000279634', gene_id)

def test_gene_id_of_invalid_name():
    with assert_raises(Exception):
        ensembl_grch38.gene_ids_of_gene_name(
            "A wonderous pony sees through your soul")

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
