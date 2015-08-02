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

from os.path import basename, join, exists
from shutil import copy2

class GenomeSource(object):
    """
    Represents a URL or local file path of a GTF or FASTA file.
    """
    def __init__(self, path_or_url, reference_name, arg_name=None):
        self.path_or_url = path_or_url
        self.reference_name = reference_name
        self.arg_name = arg_name

    @property
    def original_filename(self):
        return basename(self.path_or_url)

    @property
    def cached_filename(self):
        """
        As opposed to the remote filename or the original filename
        of a local source, this is the filename after downloading/copying
        to the cache directory.

        Unless overriden, these are the same.
        """
        return self.original_filename

    def cached_path(self, cache):
        """
        The full path to the cached file.
        """
        return join(cache.cache_directory_path,
                    self.cached_filename)

    def install_string_console(self):
        if not self.arg_name:
            raise ValueError("Expected GenomeSource to contain an arg_name: %s"
                             % str(self))

        # The shell script expects dashed arguments
        console_arg_name = self.arg_name.replace("_", "-")

        return "pyensembl install --reference-name %s --%s \"%s\"" % (
            self.reference_name, console_arg_name, self.path_or_url)

    def install_string_python(self):
        if not self.arg_name:
            raise ValueError("Expected GenomeSource to contain an arg_name: %s"
                             % str(self))
        return "Genome(reference_name=\"%s\", %s=\"%s\").install()" % (
            self.reference_name, self.arg_name, self.path_or_url)

    def is_url_format(self):
        return "://" in self.path_or_url

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

    def __str__(self):
        return "GenomeSource(path_or_url=%s, reference_name=%s, arg_name=%s)" % (
            self.path_or_url, self.reference_name, self.arg_name)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return (
            other.__class__ is GenomeSource and
            self.path_or_url == other.path_or_url and
            self.reference_name == other.reference_name and
            self.arg_name == other.arg_name)

    def __hash__(self):
        return hash((self.path_or_url, self.reference_name, self.arg_name))
