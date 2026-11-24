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

"""Implementation of the ``pyensembl verify`` CLI action."""

import json

from .formatting import format_table


def require_ensembl_releases(genomes):
    """Ensure every selected genome supports checksum verification.

    Parameters
    ----------
    genomes : list of Genome
        Genomes selected on the command line.

    Raises
    ------
    ValueError
        If any genome is neither an ``EnsemblRelease`` nor a ``PinnedGenome``.
    """
    from ..ensembl_release import EnsemblRelease
    from ..pinned_genomes import PinnedGenome

    if any(not isinstance(g, (EnsemblRelease, PinnedGenome)) for g in genomes):
        raise ValueError("verify requires an Ensembl release; use inspect for custom sources")


def run_verify(genomes, check_remote=False, json_output=False, use_color=None):
    """Verify recorded checksums of Ensembl releases without repairing them.

    Parameters
    ----------
    genomes : list of EnsemblRelease
        Releases to verify.
    check_remote : bool, optional
        Also compare recorded MD5s with current upstream manifests.
    json_output : bool, optional
        Return JSON instead of a table.
    use_color : bool, optional
        Override ANSI table header and status styling (defaults to
        :func:`pyensembl.cli.help.color_enabled` on ``sys.stdout``).

    Returns
    -------
    tuple[str, bool]
        The report text and whether every source is fine (``False`` means
        the command should exit non-zero).

    Raises
    ------
    ValueError
        If any genome is not an Ensembl or pinned release (see
        :func:`require_ensembl_releases`).
    """
    require_ensembl_releases(genomes)
    reports = [dict(release=g.release, sources=g.verify_sources(check_remote)) for g in genomes]
    if json_output:
        text = json.dumps(reports, indent=2)
    else:
        text = format_table(
            ("Release", "Source", "Local", "Upstream"),
            [
                (r["release"], s["source"], s["status"], s.get("upstream_status", "not checked"))
                for r in reports
                for s in r["sources"]
            ],
            use_color=use_color,
            status_columns=("Local", "Upstream"),
        )
    ok = not any(
        s["status"] not in ("verified", "unrecorded") or s.get("upstream_status") == "changed"
        for r in reports
        for s in r["sources"]
    )
    return text, ok
