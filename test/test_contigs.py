from pyensembl import EnsemblRelease

from test_common import releases, contigs

def test_gene_ids_on_contig():
    for release in releases:
        gene_ids_chr17 = release.gene_ids_on_contig(17)
        # gene ID of TP53
        assert "ENSG00000141510" in gene_ids_chr17, \
            (release, gene_ids_chr17[:5], len(gene_ids_chr17))

        # gene ID of SMAD4
        gene_ids_chr18 = release.gene_ids_on_contig(18)
        assert "ENSG00000141646" in gene_ids_chr18, \
            (release, gene_ids_chr18[:5], len(gene_ids_chr18))

def test_gene_names_on_contig():
    for release in releases:
        gene_names_chr17 = release.gene_names_on_contig(17)
        assert "TP53" in gene_names_chr17, \
            (release, gene_names_chr17[:5], len(gene_names_chr17))

        gene_names_chr18 = release.gene_names_on_contig(18)
        assert "SMAD4" in gene_names_chr18, \
            (release, gene_names_chr17[:5], len(gene_names_chr18))
