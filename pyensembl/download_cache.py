# Copyright (c) 2015. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from os.path import join, exists, split, abspath
from shutil import copy2, rmtree

import datacache

CACHE_BASE_SUBDIR = "pyensembl"

def cache_subdirectory(
        reference_name="",
        annotation_name="",
        annotation_version=""):
    """
    Which cache subdirectory to use for a given annotation database
    over a particular reference. All arguments can be omitted to just get
    the base subdirectory for all pyensembl cached datasets.
    """
    reference_dir = join(CACHE_BASE_SUBDIR, reference_name)
    annotation_dir = "%s%s" % (annotation_name, annotation_version)
    return join(reference_dir, annotation_dir)

class MissingRemoteFile(Exception):
    def __init__(self, url):
        self.url = url

class MissingLocalFile(Exception):
    def __init__(self, path):
        self.path = path

    def __str__(self):
        return("MissingFile(%s)" % self.path)

class DownloadCache(object):
    """
    Downloads remote files to cache, optionally copies local files into cache,
    raises custom message if data is missing.
    """
    def __init__(
            self,
            reference_name,
            annotation_name,
            annotation_version=None,
            auto_download=False,
            decompress_on_download=False,
            copy_local_files_to_cache=False,
            install_string_function=None):
        """
        Parameters
        ----------
        cache_directory_path : str
            Cache directory used as a destination for downloaded files

        reference_name : str
            Name of reference genome

        annotation_name : str
            Name of annotation database

        annotation_version : str or int
            Version or release of annotation database

        decompress_on_download : bool, optional
            If downloading a .fa.gz file, should we automatically expand it
            into a decompressed FASTA file?

        copy_local_files_to_cache : bool, optional
            If file is on the local file system, should we still copy it
            into the cache?

        cached_filename_function : fn, optional
            Function which takes the remote filename and returns a local
            renamed version

        install_string_function : fn, optional
            Function which takes dictionary of missing sources (and their
            URLs) and returns an error message with install instructions.
        """

        self.reference_name = reference_name
        self.annotation_name = annotation_name
        self.annotation_version = annotation_version

        self.cache_subdirectory = cache_subdirectory(
            reference_name=reference_name,
            annotation_name=annotation_name,
            annotation_version=annotation_version)

        # hidden since access to this variable is combined
        # with ensuring that the directpry actually exists
        self._cache_directory_path = datacache.get_data_dir(
            subdir=self.cache_subdirectory)

        self.decompress_on_download = decompress_on_download
        self.copy_local_files_to_cache = copy_local_files_to_cache

        self.install_string_function = install_string_function

    @property
    def cache_directory_path(self):
        datacache.ensure_dir(self._cache_directory_path)
        return self._cache_directory_path

    def _fields(self):
        """
        Fields used for hashing, string representation, equality comparison
        """
        return (
            ('reference_name', self.reference_name,),
            ('annotation_name', self.annotation_name),
            ('annotation_version', self.annotation_version),
            ('cache_directory_path', self.cache_directory_path),
            ('decompress_on_download', self.decompress_on_download)
            ('copy_local_files_to_cache', self.copy_local_files_to_cache)
        )

    def __eq__(self, other):
        return (
            other.__class__ is DownloadCache and
            self._fields() == other._fields()
        )

    def __hash__(self):
        return hash(self._fields())

    def __str__(self):
        fields_str = ", ".join("%s=%s" % (k, v) for (k, v) in self._fields())
        return "DownloadCache(%s)" % fields_str

    def __repr__(self):
        return str(self)

    def is_url_format(self, path_or_url):
        return "://" in path_or_url

    def cached_path(self, path_or_url):
        """
        When downloading remote files, the default behavior is to name local
        files the same as their remote counterparts.
        """
        cached_filename = split(path_or_url)[1]
        if len(cached_filename) == 0:
            raise ValueError("Can't determine local filename for %s" % (
                path_or_url,))
        return join(self.cache_directory_path, cached_filename)

    def download_or_copy_if_necessary(
            self,
            path_or_url,
            auto_download=False,
            overwrite=False):
        """
        Download a remote file or copy
        Get the local path to a possibly remote file.

        Download if file is missing from the cache directory and `auto_download`
        is True. Download even if local file exists if both `auto_download` and
        `overwrite` are True.

        If the file is on the local file system then return its path, unless
        self.copy_local_to_cache is True, and then copy it to the cache first.

        Parameters
        ----------
        path_or_url : str

        auto_download : bool, optional
            Download files if missing from local cache

        overwrite : bool, optional
            Overwrite existing copy if it exists
        """
        cached_path = self.cached_path(path_or_url)

        if exists(cached_path) and not overwrite:
            return cached_path

        if self.is_url_format(path_or_url):
            if auto_download:
                cached_filename = split(cached_path)[1]
                return datacache.fetch_file(
                    download_url=path_or_url,
                    filename=cached_filename,
                    decompress=self.decompress_on_download,
                    subdir=self.cache_subdirectory,
                    force=overwrite)
            else:
                raise MissingRemoteFile(path_or_url)
        else:
            local_path = abspath(path_or_url)
            if not exists(local_path):
                raise MissingLocalFile(local_path)
            elif self.copy_local_files_to_cache:
                copy2(local_path, cached_path)
                return cached_path
            else:
                return local_path

    def _raise_missing_file_error(self, missing_urls_dict):
        if self.install_string_function is None:
            return "Missing genome data files from %s" % missing_urls_dict

        install_string = self.install_string_function(missing_urls_dict)
        missing_urls = list(missing_urls_dict.values())
        assert len(missing_urls_dict) > 0
        if len(missing_urls) == 1:
            raise ValueError("Missing genome data file from %s, run: %s" % (
                missing_urls[0], install_string))
        else:
            raise ValueError("Missing genome data files from %s, run: %s" % (
                missing_urls, install_string))

    def local_path_or_install_error(
            self,
            field_name,
            path_or_url,
            auto_download=False,
            overwrite=False):
        try:
            return self.download_or_copy_if_necessary(
                path_or_url,
                auto_download=auto_download,
                overwrite=overwrite)
        except MissingRemoteFile:
            self._raise_missing_file_error({field_name: path_or_url})

    def local_paths_or_install_error(
            self,
            sources_dict,
            auto_download=False,
            overwrite=False):
        """
        Constructs result dictionary with local path for each (k, path_or_url)
        mapping in the sources dict.
        """
        missing_urls_dict = {}
        results = {}
        for (field_name, path_or_url) in sources_dict.items():
            try:
                results[field_name] = self.local_path_or_install_error(
                    path_or_url)
            except MissingRemoteFile:
                missing_urls_dict[field_name] = path_or_url

        if len(missing_urls_dict) == 0:
            return results

        self._raise_missing_file_error(missing_urls_dict)

    def delete_all_files(self):
        rmtree(self.cache_directory_path)
