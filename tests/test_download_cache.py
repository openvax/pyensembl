import gzip
import os
import stat
import tempfile

import pytest
from pytest import raises as assert_raises
from pyensembl.download_cache import DownloadCache, MissingLocalFile, MissingRemoteFile


from .data import data_path
from .common import ok_

download_cache = DownloadCache(
    reference_name="__test_reference",
    annotation_name="__test_annotation",
    copy_local_files_to_cache=False,
)


def test_download_cache_missing_local_file():
    # clear the cache
    download_cache.delete_cache_directory()
    with assert_raises(MissingLocalFile):
        download_cache.download_or_copy_if_necessary(
            path_or_url="test_file_doesn_not_exist.file"
        )


def test_download_cache_missing_remote_file():
    # clear the cache
    download_cache.delete_cache_directory()
    with assert_raises(MissingRemoteFile):
        download_cache.download_or_copy_if_necessary(
            path_or_url="ftp://NOTAURL.NOTAURL.NOTAURL"
        )


def test_download_cache_custom_location():
    test_file = "refseq.ucsc.small.gtf"
    tmp_dir = tempfile.gettempdir()

    print("DIR: %s" % tmp_dir)
    assert tmp_dir is not None

    os.environ["PYENSEMBL_CACHE_DIR"] = tmp_dir

    # We need another instance of DownloadCache
    # that copies files over to cache folder
    download_cache = DownloadCache(
        reference_name="test_reference",
        annotation_name="test_annotation",
        copy_local_files_to_cache=True,
    )

    # clean up
    download_cache.delete_cache_directory()
    download_cache.download_or_copy_if_necessary(
        download_if_missing=True, path_or_url=data_path(test_file)
    )

    full_path = os.path.join(
        tmp_dir, "pyensembl", "test_reference", "test_annotation", test_file
    )
    print("FULL PATH: %s" % full_path)
    assert len(full_path) > 0

    ok_(os.path.exists(full_path))
    del os.environ["PYENSEMBL_CACHE_DIR"]


@pytest.mark.parametrize("decompress", [False, True])
def test_downloads_follow_the_umask_and_decompress_on_request(tmp_path, decompress):
    source = tmp_path / "annotation.gtf.gz"
    source.write_bytes(gzip.compress(b"gtf rows\n"))
    cache = DownloadCache(
        "reference", "annotation",
        decompress_on_download=decompress,
        cache_directory_path=str(tmp_path / "cache"),
    )
    previous = os.umask(0o002)
    try:
        path = cache.download_or_copy_if_necessary(source.as_uri(), download_if_missing=True)
    finally:
        os.umask(previous)
    # Shared caches must be readable by the group, not only the downloader.
    assert stat.S_IMODE(os.stat(path).st_mode) == 0o664
    assert path.endswith(".gtf") == decompress
    with open(path, "rb") as handle:
        data = handle.read()
    assert (data if decompress else gzip.decompress(data)) == b"gtf rows\n"
