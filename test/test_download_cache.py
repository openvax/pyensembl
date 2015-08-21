from __future__ import absolute_import

from nose.tools import assert_raises
from pyensembl.download_cache import (
    DownloadCache,
    MissingLocalFile,
    MissingRemoteFile
)

download_cache = DownloadCache(
    reference_name="__test_reference",
    annotation_name="__test_annotation",
    copy_local_files_to_cache=False)

def test_download_cache_missing_local_file():
    # clear the cache
    download_cache.delete_cache_directory()
    with assert_raises(MissingLocalFile):
        download_cache.download_or_copy_if_necessary(
            path_or_url="test_file_doesn_not_exist.file")

def test_download_cache_missing_remote_file():
    # clear the cache
    download_cache.delete_cache_directory()
    with assert_raises(MissingRemoteFile):
        download_cache.download_or_copy_if_necessary(
            path_or_url="ftp://NOTAURL.NOTAURL.NOTAURL")
