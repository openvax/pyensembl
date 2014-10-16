from pyensembl import EnsemblRelease

ensembl = EnsemblRelease(77)
# chr6:29,945,884  is a position for HLA-A
# based on:
# http://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000206503;r=6:29941260-29945884

def test_gene_name_hla_a():
    # chr6:29,945,884  is a position for HLA-A
    # based on:
    # http://useast.ensembl.org/Homo_sapiens/Gene/
    # Summary?db=core;g=ENSG00000206503;r=6:29941260-29945884
    names = ensembl.gene_names_at_locus(6, 29945884)
    assert names == ["HLA-A"], names

def test_gene_ids_hla_a():
    ids = ensembl.gene_ids_at_locus(6, 29945884)
    assert ids == ['ENSG00000206503'], ids
