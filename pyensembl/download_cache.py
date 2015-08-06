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
from shutil import copy2

import datacache

CACHE_BASE_SUBDIR = "pyensembl"

def cache_subdirectory(
        annotation_name=None,
        annotation_version=None,
        reference_name=None):
    result = CACHE_BASE_SUBDIR
    for subdir in [annotation_name, annotation_version, reference_name]:
        if subdir is not None:
            result = join(result, str(subdir))
    return result

class MissingRemoteFile(Exception):
    def __init__(self, path_or_url):
        self.path_or_url = path_or_url

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
            force_download=False,
            decompress_on_download=False,
            copy_local_to_cache=False,
            local_filename_function=None,
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

        auto_download : bool, optional
            Download files if missing from local cache

        force_download : bool, optional
            Download files even if already cached

        decompress_on_download : bool, optional
            If downloading a .fa.gz file, should we automatically expand it
            into a decompressed FASTA file?

        copy_local_to_cache : bool, optional
            If file is on the local file system, should we still copy it
            into the cache?

        local_filename_function : fn, optional
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

        self.auto_download = auto_download
        self.force_download = force_download
        self.decompress_on_download = decompress_on_download
        self.copy_local_to_cache = copy_local_to_cache

        if local_filename_function:
            self.local_filename_function = local_filename_function
        else:
            self.local_filename_function = self.default_local_filename_function

        if install_string_function:
            self.install_string_function = install_string_function
        else:
            self.install_string_function = self.default_install_string

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
            ('auto_download', self.auto_download),
            ('force_download', self.force_download),
            ('copy_local_to_cache', self.copy_local_to_cache),
            ('decompress_on_download', self.decompress_on_download)
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

    def default_local_filename_function(self, path_or_url):
        """
        Default is to name local files the same as their remote counterparts.
        """
        remote_filename = split(path_or_url)[1]
        if len(remote_filename) == 0:
            raise ValueError("Can't determine local filename for %s" % (
                path_or_url,))
        return remote_filename

    def default_install_string(self, missing_files):
        """
        Add every missing file to the install string shown to the user
        in an error message.
        """
        args = [
            "--reference_name", self.reference_name,
            "--annotation_name", self.annotation_name]
        if self.annotation_version:
            args.extend(["--annotation-version", str(self.annotation_version)])
        for (name, url) in missing_files:
            args.append("--%s" % (name.replace("_", "-"),))
            args.append("\"%s\"" % (url,))

        return "pyensembl install %s" % " ".join(args)

    def is_url_format(self, path_or_url):
        return "://" in path_or_url

    def _local_path(self, path_or_url):
        """
        Get the local path to a possibly remote file.

        Always download remote  files if self.force_download is True or if the
        file is missing from the cache directory and self.auto_download is True.
        If the file is on the local file system then return its path, unless
        self.copy_local_to_cache is True, and then copy it to the cache first.
        """
        cache_filename = self.local_filename_function(path_or_url)
        cached_path = join(self.cache_directory_path, cache_filename)

        if exists(cached_path) and not self.force_download:
            return cached_path
        if self.is_url_format(path_or_url):
            if self.auto_download or self.force_download:
                return datacache.fetch_file(
                    download_url=path_or_url,
                    filename=cache_filename,
                    decompress=self.decompress_on_download,
                    subdir=self.cache_subdirectory,
                    force=self.force_download)
            else:
                raise MissingRemoteFile(path_or_url)
        else:
            local_path = abspath(path_or_url)
            if not exists(local_path):
                raise ValueError(
                    "Couldn't find genome data file %s" % local_path)
            elif self.copy_local_to_cache:
                copy2(local_path, cached_path)
                return cached_path
            else:
                return local_path
        return local_path

    def raise_missing_file_error(self, missing_urls_dict):
        install_string = self.install_string_function(missing_urls_dict)
        missing_urls = list(missing_urls_dict.values())
        assert len(missing_urls_dict) > 0
        if len(missing_urls) == 1:
            raise ValueError("Missing genome data file from %s, run: %s" % (
                missing_urls[0], install_string))
        else:
            raise ValueError("Missing genome data files from %s, run: %s" % (
                missing_urls, install_string))

    def local_path(self, field_name, path_or_url):
        try:
            return self._local_path(path_or_url)
        except MissingRemoteFile:
            self.raise_missing_file_error({field_name: path_or_url})

    def local_paths(self, **sources):
        """
        Constructs result dictionary with local path for each (k, path_or_url)
        mapping in the sources dict.
        """
        missing_urls_dict = {}
        results = {}
        for (field_name, path_or_url) in sources.items():
            try:
                local_path = self._local_path(path_or_url)
                results[field_name] = local_path
            except MissingRemoteFile:
                missing_urls_dict[field_name] = path_or_url

        if len(missing_urls_dict) == 0:
            return results

        self.raise_missing_file_error(missing_urls_dict)