from __future__ import absolute_import

from .common import test_ensembl_releases
from .data import TP53_gene_id

@test_ensembl_releases()
def test_TP53_gene_object_by_id(ensembl):
    # when we look up TP53 by its gene ID, we should get the
    # correct gene back
    gene = ensembl.gene_by_id(TP53_gene_id)
    assert gene.name == "TP53", \
        "Incorrect gene name %s for gene ID %s in %s" % (
            gene.name, gene.id, ensembl)
    assert gene.contig == "17", \
        "Incorrect gene contig %s for gene ID %s in %s" % (
            gene.contig, gene.id, ensembl)

@test_ensembl_releases()
def test_TP53_gene_object_by_name(ensembl):
    genes = ensembl.genes_by_name("TP53")
    # we should only have one TP53 gene (there aren't any copies)
    assert len(genes) == 1, \
        "Expected only one gene with name TP53, got %s" % (genes,)
    # make sure it has the correct gene ID
    assert genes[0].id == TP53_gene_id, \
        "Expected gene to have ID %s, got %s" % (TP53_gene_id, genes[0].id)

@test_ensembl_releases()
def test_equal_genes(ensembl):
    gene1 = ensembl.genes_by_name("TP53")[0]
    # get an identical gene
    gene2 = ensembl.gene_by_id(gene1.id)

    assert hash(gene1) == hash(gene2)
    assert gene1 == gene2

@test_ensembl_releases()
def test_not_equal_genes(release):
    gene1 = release.genes_by_name("MUC1")[0]
    gene2 = release.genes_by_name("BRCA1")[0]
    assert hash(gene1) != hash(gene2)
    assert gene1 != gene2
