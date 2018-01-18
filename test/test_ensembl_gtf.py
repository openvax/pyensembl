from __future__ import absolute_import
from os.path import exists


from .common import test_ensembl_releases

@test_ensembl_releases()
def gtf_path_endswith_gtf_gz(ensembl):
    path = ensembl.gtf.gtf_path
    assert exists(path)
    assert path.endswith(".gtf.gz")
