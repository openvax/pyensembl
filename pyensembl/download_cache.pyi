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

from typing import TYPE_CHECKING, Dict, List, Optional, Tuple, Union

if TYPE_CHECKING:
    import logging

logger: logging.Logger = ...

CACHE_BASE_SUBDIR = "pyensembl"
CACHE_DIR_ENV_KEY = "PYENSEMBL_CACHE_DIR"

def cache_subdirectory(
    reference_name: Optional[str] = None,
    annotation_name: Optional[str] = None,
    annotation_version: Optional[Union[str, int]] = None,
) -> str: ...

class MissingRemoteFile(Exception):
    def __init__(self, url: str) -> None: ...

class MissingLocalFile(Exception):
    def __init__(self, path: str) -> None: ...
    def __str__(self) -> str: ...

class DownloadCache(object):
    def __init__(
        self,
        reference_name: str,
        annotation_name: str,
        annotation_version: Union[str, int] = None,
        decompress_on_download: bool = False,
        copy_local_files_to_cache: bool = False,
        install_string_function: Optional[function] = None,
        cache_directory_path: Optional[str] = None,
    ) -> None: ...
    @property
    def cache_directory_path(self) -> str: ...
    def _fields(self) -> Tuple[Tuple[str, Union[str, int, bool]]]: ...
    def __eq__(self, other) -> bool: ...
    def __hash__(self) -> int: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
    def is_url_format(self, path_or_url: str) -> bool: ...
    def _remove_compression_suffix_if_present(self, filename: str) -> str: ...
    def cached_path(self, path_or_url: str) -> str: ...
    def _download_if_necessary(
        self, url: str, download_if_missing: bool, overwrite: bool
    ) -> str: ...
    def _copy_if_necessary(self, local_path: str, overwrite: bool) -> str: ...
    def download_or_copy_if_necessary(
        self,
        path_or_url: str,
        download_if_missing: bool = False,
        overwrite: bool = False,
    ) -> str: ...
    def _raise_missing_file_error(self, missing_urls_dict: Dict) -> None: ...
    def local_path_or_install_error(
        self,
        field_name: str,
        path_or_url: str,
        download_if_missing: bool = False,
        overwrite: bool = False,
    ) -> str: ...
    def delete_cached_files(
        self, prefixes: List[str] = [], suffixes: List[str] = []
    ) -> None: ...
    def delete_cache_directory(self) -> None: ...
