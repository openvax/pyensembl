from test_common import releases


def test_TP53_gene_object_by_id():
    for release in releases:
        print release
        # when we look up TP53 by its gene ID, we should get the
        # correct gene back
        gene = release.gene_by_id("ENSG00000141510")
        assert gene.name == "TP53", gene
        assert gene.contig == "17", gene

def test_TP53_gene_object_by_name():
    for release in releases:
        print release
        genes = release.genes_by_name("TP53")
        # we should only have one TP53 gene (there aren't any copies)
        assert len(genes) == 1
        # make sure it has the correct gene ID
        assert genes[0].id == "ENSG00000141510", genes
