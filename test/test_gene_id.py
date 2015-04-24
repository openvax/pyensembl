from __future__ import absolute_import

from pyensembl import ensembl77 as ensembl

from nose.tools import assert_raises, ok_

def test_gene_ids_of_gene_name_hla_release77():
    hla_a_gene_ids = ensembl.gene_ids_of_gene_name("HLA-A")
    assert 'ENSG00000206503' in hla_a_gene_ids, hla_a_gene_ids

    hla_b_gene_ids = ensembl.gene_ids_of_gene_name("HLA-B")
    assert 'ENSG00000234745' in hla_b_gene_ids, hla_b_gene_ids

    hla_c_gene_ids = ensembl.gene_ids_of_gene_name("HLA-C")
    assert 'ENSG00000204525' in hla_c_gene_ids, hla_c_gene_ids

def test_gene_id_of_protein_id_release77():
    gene_id = ensembl.gene_id_of_protein_id("ENSP00000485677")
    ok_('ENSG00000279634', gene_id)

def test_gene_id_of_invalid_name():
    with assert_raises(Exception, None):
        ensembl.gene_ids_of_gene_name("A wonderous pony sees through your soul")
