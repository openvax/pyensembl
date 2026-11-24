# PYTHON_ARGCOMPLETE_OK
# (argcomplete's global activation completes only scripts with this marker.)
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

"""Compatibility entry point for the ``pyensembl`` console script.

The CLI implementation lives in :mod:`pyensembl.cli` (organized by concern:
``main``, ``parser``, ``formatting``, ``logging_config``, ``selection``,
``install``, ``listing``, ``available``, ``inspect``, ``verify``, and
``delete``). This module keeps ``pyensembl.shell:run`` and the ``pyensembl``
console script importable while keeping module-level imports cheap for
``argcomplete`` tab completion.
"""

from .cli.available import (  # noqa: F401
    available_catalogue,
    format_available_species,
    run_available,
    selects_catalogue,
)
from .cli.delete import delete_genome_files, run_prune  # noqa: F401
from .cli.inspect import inspect_genomes  # noqa: F401
from .cli.install import (  # noqa: F401
    InstallError,
    install_genome,
    progress_available,
    run_install,
)
from .cli.listing import (  # noqa: F401
    collect_all_installed_ensembl_releases,
    format_cache_dirs,
    format_installed_genomes,
    format_installed_genomes_json,
    genome_fasta_status,
)
from .cli.logging_config import configure_logging  # noqa: F401
from .cli.main import main, run  # noqa: F401
from .cli.parser import (  # noqa: F401
    DELETE_ACTIONS,
    build_parser,
    complete_release,
    complete_species,
    installed_releases,
    parser,
    validate_arguments,
)
from .cli.selection import (  # noqa: F401
    all_combinations_of_ensembl_genomes,
    collect_selected_genomes,
)
from .cli.verify import require_ensembl_releases, run_verify  # noqa: F401
