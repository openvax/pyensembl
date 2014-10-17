
from pyensembl import EnsemblRelease

releases = [
    EnsemblRelease(75),
    EnsemblRelease(77)
]

contigs = range(1,23) + ["X", "Y", "M"]

def test_gene_id_length_on_contig():
    for release in releases:
        for contig in contigs:
            gene_ids = release.gene_ids_on_contig(contig)
            # Ensembl gene ids are formatted like ENSG00000223972
            # which is always length 15
            for gene_id in gene_ids:
                assert len(gene_id) == 15, \
                    "Gene ID with wrong length: %s on contig %s of %s" % (
                        gene_id,
                        contig,
                        release
                    )