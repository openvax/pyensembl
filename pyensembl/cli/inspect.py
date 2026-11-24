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

"""Implementation of the ``pyensembl inspect`` CLI action."""

from dataclasses import fields
import json

from .formatting import display_path, format_table, genome_description


def inspect_genomes(genomes, check_genome_fasta=False, json_output=False, use_color=None):
    """Render the offline file inventory exposed by ``Genome.inspect_data``.

    Parameters
    ----------
    genomes : sequence of Genome
        Selected genomes to inspect offline.
    check_genome_fasta : bool, optional
        When ``True``, verify reference DNA index integrity on disk.
    json_output : bool, optional
        When ``True``, return formatted JSON instead of human-readable tables.
    use_color : bool, optional
        Override ANSI table header and status styling (defaults to
        :func:`pyensembl.cli.help.color_enabled` on ``sys.stdout``).

    Returns
    -------
    str
        Formatted table(s) or JSON array describing each genome's files.
    """
    reports = []
    tables = []
    for genome in genomes:
        report = genome.inspect_data(check_genome_fasta=check_genome_fasta)
        report.update(
            reference_name=genome.reference_name,
            annotation_name=genome.annotation_name,
            annotation_version=genome.annotation_version,
            cache_directory=genome.download_cache.cache_directory_path,
        )
        if json_output:
            report["files"] = {
                role: {
                    field.name: (str(file.error) if file.error else None)
                    if field.name == "error"
                    else getattr(file, field.name)
                    for field in fields(file)
                }
                for role, file in report["files"].items()
            }
            reports.append(report)
        else:
            rows = [
                (
                    role,
                    file.status,
                    str(file.size) if file.size is not None else "-",
                    display_path(file.path),
                    file.source_url or "-",
                )
                for role, file in report["files"].items()
            ]
            lines = [
                "%s: annotation=%s, reference DNA=%s"
                % (
                    genome_description(genome),
                    report["annotation"] or "not configured",
                    report["reference_dna"] or "not configured",
                ),
                format_table(
                    ("File", "Status", "Bytes", "Path", "Fetched from"),
                    rows,
                    use_color=use_color,
                    status_columns=("Status", "Bytes", "Fetched from"),
                ),
            ]
            lines.extend(
                "%s: %s" % (role, file.error)
                for role, file in report["files"].items()
                if file.error
            )
            tables.append("\n".join(lines))
    return json.dumps(reports, indent=2) if json_output else "\n\n".join(tables)
