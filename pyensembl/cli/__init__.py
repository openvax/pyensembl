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

"""Command-line interface package for PyEnsembl.

Imports stay cheap at package load time so tab completion via ``argcomplete``
does not import genome, database, or download modules.
"""

from .help import (
    HelpArgumentParser,
    HelpFormatter,
    color_enabled,
    style_help,
    style_status,
    style_table_header,
)
from .main import run
from .parser import build_parser

__all__ = [
    "HelpArgumentParser",
    "HelpFormatter",
    "build_parser",
    "color_enabled",
    "run",
    "style_help",
    "style_status",
    "style_table_header",
]
