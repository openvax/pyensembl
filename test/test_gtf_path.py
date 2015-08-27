from __future__ import absolute_import

from glob import glob
from os.path import join
from pyensembl import GTF

from .data import data_path
from .common import TemporaryDirectory

gtf_path = data_path("mouse.ensembl.81.partial.ENSMUSG00000017167.gtf")

def test_gtf_object_path():
    """
    Make sure that GTF doesn't do anything funny under the hood
    (such as copying local files to some cached directory), all of that
    should be done at the level of Genome/DownloadCache.
    """
    gtf_object = GTF(gtf_path)
    assert gtf_path == gtf_object.gtf_path

def test_gtf_creates_csv_files_in_cache_dir():
    """
    Make sure that GTF obeys the cache_dir argument by creating all of its
    index files there.
    """
    with TemporaryDirectory() as tmpdir:
        gtf_object = GTF(gtf_path, cache_directory_path=tmpdir)
        csv_path = gtf_object.data_subset_path()
        assert tmpdir in csv_path, \
            "Expected CSV path %s to contain cache_dir" % csv_path
        search_pattern = join(tmpdir, "*mouse*")
        assert len(glob(search_pattern)) == 0, \
            "Temporary directory should be empty"
        # creating the dataframe should have the effect of triggering
        # GTF parsing and then saving the parsed results in a csv file
        gtf_object.dataframe()
        assert len(glob(search_pattern)) > 0, \
            "Expected GTF to save files in cache_dir"