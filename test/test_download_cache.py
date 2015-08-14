from __future__ import absolute_import

from pyensembl.download_cache import DownloadCache, MissingLocalFile

download_cache = DownloadCache(
    reference_name="__test_reference",
    annotation_name="__test_annotation",
    copy_local_files_to_cache=True)

def test_download_cache_missing_file():
    # clear the cache
    download_cache.delete_all_files()
    try:
        download_cache.download_or_copy_if_necessary(
            path_or_url="test_file_doesn_not_exist.file")
        assert False, "Should not succeed"
    except MissingLocalFile:
        pass
