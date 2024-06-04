# Copyright (c) 2015-2016. Mount Sinai School of Medicine
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

from typing import TYPE_CHECKING, Dict, Generator, Tuple, Union

if TYPE_CHECKING:
    import logging
    from io import BufferedIOBase

logger: logging.Logger = ...

def _parse_header_id(line: bytes) -> str: ...

class FastaParser(object):
    def __init__(self) -> None: ...
    def read_file(self, fasta_path: str) -> Dict[str, str]: ...
    def iterate_over_file(
        self, fasta_path: str
    ) -> Generator[Tuple[str, str], None, None]: ...
    def _open(self, fasta_path: str) -> Union[BufferedIOBase]: ...
    def _current_entry(self) -> Tuple[str, str]: ...
    def _read_header(self, line: bytes) -> Tuple[str, str]: ...

def parse_fasta_dictionary(fasta_path: str) -> Dict[str, str]: ...
