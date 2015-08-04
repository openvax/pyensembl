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

from __future__ import print_function, absolute_import

from os.path import join
import datacache

CACHE_BASE_SUBDIR = "pyensembl"

def cache_subdir(annotation_name, annotation_version=None, reference_name=None):
    result = CACHE_BASE_SUBDIR
    for subdir in [annotation_name, annotation_version, reference_name]:
        if subdir is not None:
            result = join(result, str(subdir))
    return result

def get_download_cache(
        annotation_name,
        annotation_version=None,
        reference_name=None):
    subdir = cache_subdir(annotation_name, annotation_version, reference_name)
    return datacache.Cache(subdir)
