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

from os.path import basename, join, exists, split, abspath
from shutil import copy2

import datacache

CACHE_BASE_SUBDIR = "pyensembl"

def cache_subdir(
        annotation_name=None,
        annotation_version=None,
        reference_name=None):
    result = CACHE_BASE_SUBDIR
    for subdir in [annotation_name, annotation_version, reference_name]:
        if subdir is not None:
            result = join(result, str(subdir))
    return result

def cache_directory_path(
        annotation_name=None,
        annotation_version=None,
        reference_name=None):
    """
    What directory are cached files placed in for a particular genome
    annotation?
    """
    subdir = cache_subdir(
        annotation_name=annotation_name,
        annotation_version=annotation_version,
        reference_name=reference_name)
    return datacache.data_dir(subdir=subdir)

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

        self.cache_directory_path = cache_directory_path(
            reference_name=reference_name,
            annotation_name=annotation_name,
            annotation_version=annotation_version)

        self.auto_download = auto_download
        self.force_download = force_download
        self.copy_local_to_cache = copy_local_to_cache
        if local_filename_function:
            self.local_filename_function = local_filename_function
        else:
            self.local_filename_function = self.default_local_filename_function
        if install_string_function:
            self.install_string_function = install_string_function
        else:
            self.install_string_function = self.default_install_string

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
            ('copy_local_to_cache', self.copy_local_to_cache)
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

    def default_local_filename_function(self, url):
        """
        Default is to name local files the same as their remote counterparts.
        """
        remote_filename = split(url)[1]
        if len(remote_filename) == 0:
            raise ValueError("Can't determine local filename for %s" % (
                url,))
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

    def local_path(self, path_or_url):
        if self.is_url(path_or_url):
            local_filename = self.local_filename_function(path_or_url)
            local_path = join(self.cache_directory_path, local_filename)
        else:
            if self.copy_local_to_cache:
                # COPY
                local_path
        if not exists(local_path):
            raise ValueError("Missing local file: %s" % (local_path,))
        return local_path

    def local_paths(self, **sources):
        """
        Constructs result dictionary with local path for each (k, path_or_url)
        mapping in the sources dict.
        """
        return {
            (k, self.local_path(path_or_url))
            for (k, path_or_url) in sources.items()
        }

    def copy_to_cache_if_needed(self, cache, force):
        """
        Given a path to a file, copy the file to the cache's directory
        and return the path to the new file within the cache's directory.

        If `force` is True, overwrites an existing cache directory file.
        """
        assert not self.is_url_format(), \
            "Copying should not be called on a URL: %s" % self.path_or_url
        new_path = self.cached_path(cache)
        if not exists(new_path) or force:
            copy2(self.path_or_url, new_path)
        return new_path
