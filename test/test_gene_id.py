from pyensembl import EnsemblRelease

from nose.tools import raises

ensembl = EnsemblRelease(77)

def test_gene_id_of_gene_name_hla_release77():
    assert ensembl.gene_id_of_gene_name("HLA-A") == 'ENSG00000206503'
    assert ensembl.gene_id_of_gene_name("HLA-B") == 'ENSG00000234745'
    assert ensembl.gene_id_of_gene_name("HLA-C") == 'ENSG00000204525'

@raises(Exception)
def test_gene_id_of_invalid_name():
    ensembl.gene_id_of_gene_name("A wonderous pony sees through your soul")




