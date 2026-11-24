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

"""Implementation of the ``pyensembl available`` CLI action."""

import json
import sys

from .formatting import (
    DIVISION_LABELS,
    format_release_range,
    format_table,
    species_display_name,
)
from .help import (
    ANSI_BOLD_CYAN,
    ANSI_DIM,
    ANSI_RESET,
    color_enabled,
    style_status,
)


def selects_catalogue(
    *,
    species=(),
    assembly=None,
    assembly_accession=None,
    provider=None,
    reference_name=None,
    include_alt=None,
    json_output=False,
):
    """Return whether ``available`` should render the filtered release catalogue.

    Parameters
    ----------
    species : sequence of str, optional
        Species filter values from ``--species``.
    assembly : str, optional
        Assembly filter from ``--assembly``.
    assembly_accession : str, optional
        Versioned accession filter from ``--assembly-accession``.
    provider : str, optional
        Dated-release provider filter from ``--provider``.
    reference_name : str, optional
        Reference name selector from ``--reference-name``.
    include_alt : bool, optional
        Explicit dated annotation coverage filter from ``--include-alt``.
    json_output : bool, optional
        Whether ``--json`` was requested.

    Returns
    -------
    bool
        ``True`` when any selector flag or ``json_output`` is active,
        ``False`` for the default table of supported species.
    """
    return bool(
        species
        or assembly
        or assembly_accession
        or provider
        or reference_name
        or include_alt is not None
        or json_output
    )


def available_catalogue(
    *,
    species=(),
    assembly=None,
    assembly_accession=None,
    provider=None,
    reference_name=None,
    include_alt=None,
    refresh=False,
    json_output=False,
):
    """Render the filtered Ensembl/pinned release catalogue as text or JSON.

    Parameters
    ----------
    species : sequence of str, optional
        Species filter values.
    assembly : str, optional
        Assembly name/alias selector.
    assembly_accession : str, optional
        Explicit assembly accession selector.
    provider : str, optional
        Annotation provider filter.
    reference_name : str, optional
        Reference name selector (cannot conflict with ``assembly``).
    include_alt : bool, optional
        Explicit dated coverage selector.
    refresh : bool, optional
        When ``True``, refresh dated listings from Ensembl.
    json_output : bool, optional
        When ``True``, return formatted JSON instead of a table.

    Returns
    -------
    str
        Formatted table or JSON array of matching release records.

    Raises
    ------
    ValueError
        If both ``reference_name`` and ``assembly`` are given and disagree.
    """
    from ..release_catalogue import available_releases

    if reference_name and assembly and reference_name != assembly:
        raise ValueError("Use one assembly selector with available")
    rows = []
    for species_name in species or [None]:
        rows.extend(
            available_releases(
                species=species_name,
                assembly=assembly or reference_name,
                assembly_accession=assembly_accession,
                provider=provider,
                include_alt=include_alt,
                refresh=refresh,
            )
        )
    if json_output:
        return json.dumps(rows, indent=2)
    return format_table(
        ("Release", "Assembly", "Provider", "Coverage", "Evidence", "Latest"),
        [
            (
                r["release"],
                r.get("assembly_name") or r.get("assembly_accession"),
                r["provider"],
                r["coverage"],
                "subset" if r["subset"] else (
                    "; ".join(
                        c["relationship"] + " to " + str(c["numbered"])
                        for c in r["comparisons"]
                    )
                    or (
                        "frozen"
                        if r["kind"] == "numbered"
                        else ("pinned" if r["kind"] == "pinned" else "uncompared")
                    )
                ),
                "yes" if r["latest"] else "",
            )
            for r in rows
        ],
    )


def format_available_species(use_color=None, dated_releases=None):
    """Render the overview table of supported species and Ensembl release ranges.

    When ``use_color`` is ``None`` (the default), ANSI styling is applied
    according to :func:`pyensembl.cli.help.color_enabled` on ``sys.stdout``.
    ``dated_releases`` optionally maps latin names to the annotation dates of
    their registered dated assembly (``None`` when unknown) and adds a column
    for them.

    Parameters
    ----------
    use_color : bool, optional
        Override ANSI styling (defaults to
        :func:`pyensembl.cli.help.color_enabled` on ``sys.stdout``).
    dated_releases : dict of str to (list of str or None), optional
        Optional dated-release dates keyed by species latin name.

    Returns
    -------
    str
        Formatted division-grouped species table.
    """
    from ..pinned_genomes import pinned_assemblies
    from ..species import Species

    if use_color is None:
        use_color = color_enabled(sys.stdout)
    heading = ANSI_BOLD_CYAN if use_color else ""
    dim = ANSI_DIM if use_color else ""
    reset = ANSI_RESET if use_color else ""

    species_by_division = {key: [] for key, _ in DIVISION_LABELS}
    for species in sorted(Species.all(), key=lambda s: s.latin_name):
        species_by_division.setdefault(species.division, []).append(species)

    all_species = [s for group in species_by_division.values() for s in group]
    if not all_species:
        return ""

    def column_width(values, fallback):
        return max((len(v) for v in values), default=fallback)

    pinned_items = pinned_assemblies()
    pinned_by_species = {}
    for item in pinned_items:
        pinned_by_species.setdefault(item.species, []).append(item)

    name_w = column_width([species_display_name(s) for s in all_species], 8)
    asm_w = column_width(
        [asm for s in all_species for asm in s.reference_assemblies]
        + [item.reference_name for item in pinned_items],
        8,
    )
    rng_w = column_width(
        [
            format_release_range(start, end)
            for s in all_species
            for (start, end) in s.reference_assemblies.values()
        ] + ["%s (pinned)" % item.annotation_version for item in pinned_items],
        8,
    )
    latin_w = column_width([s.latin_name for s in all_species], 8)

    def dates_cell(species, assembly):
        if dated_releases is None or species.dated_releases is None:
            return ""
        if assembly != species.dated_reference_name:
            return ""
        dates = dated_releases.get(species.latin_name)
        return "?" if dates is None else ", ".join(dates)

    dates_w = column_width(
        [dates_cell(s, asm) for s in all_species for asm in s.reference_assemblies], 0
    )

    col_name = max(name_w, len("Species")) + 2
    col_asm = max(asm_w, len("Assembly")) + 2
    col_rng = max(rng_w, len("Releases")) + 2
    col_dates = max(dates_w, len("Dated releases")) + 2 if dated_releases is not None else 0
    col_latin = max(latin_w, len("Latin name"))
    total_w = col_name + col_asm + col_rng + col_dates + col_latin

    lines = []
    header_row = "%-*s%-*s%-*s%-*s%s" % (
        col_name, "Species",
        col_asm, "Assembly",
        col_rng, "Releases",
        col_dates, "Dated releases" if col_dates else "",
        "Latin name",
    )
    lines.append("%s%s%s" % (heading, header_row, reset))
    lines.append("─" * total_w)

    for division_key, division_label in DIVISION_LABELS:
        members = species_by_division.get(division_key, [])
        if not members:
            continue
        members.sort(key=species_display_name)
        lines.append("")
        section = "── %s " % division_label
        section += "─" * max(2, total_w - len(section))
        lines.append("%s%s%s" % (heading, section, reset))
        for species in members:
            for i, (asm, (start, end)) in enumerate(
                species.reference_assemblies.items()
            ):
                if i == 0:
                    name_cell = species_display_name(species)
                    latin_cell = "%s%s%s" % (dim, species.latin_name, reset)
                else:
                    name_cell = ""
                    latin_cell = ""
                dates_raw = dates_cell(species, asm)
                dates_rendered = (
                    style_status(dates_raw, use_color=use_color)
                    + (" " * (col_dates - len(dates_raw)))
                    if col_dates
                    else ""
                )
                lines.append(
                    "%-*s%-*s%-*s%s%s"
                    % (
                        col_name, name_cell,
                        col_asm, asm,
                        col_rng, format_release_range(start, end),
                        dates_rendered,
                        latin_cell,
                    )
                )
            for item in pinned_by_species.get(species.latin_name, ()):
                lines.append(
                    "%-*s%-*s%-*s%-*s%s"
                    % (
                        col_name, "",
                        col_asm, item.reference_name,
                        col_rng, "%s (pinned)" % item.annotation_version,
                        col_dates, "",
                        "",
                    )
                )
    if any(
        dates_cell(species, asm) == "?"
        for species in all_species for asm in species.reference_assemblies
    ):
        lines.append("")
        lines.append("? = dates could not be fetched and none are cached")
    return "\n".join(lines)


def run_available(
    *,
    species=(),
    assembly=None,
    assembly_accession=None,
    provider=None,
    reference_name=None,
    include_alt=None,
    refresh=False,
    json_output=False,
):
    """Describe available Ensembl species, assemblies, and releases.

    Parameters
    ----------
    species : sequence of str, optional
        Species filter values.
    assembly : str, optional
        Assembly name/alias selector.
    assembly_accession : str, optional
        Versioned assembly accession selector.
    provider : str, optional
        Annotation provider filter.
    reference_name : str, optional
        Reference name selector.
    include_alt : bool, optional
        Explicit dated annotation coverage filter.
    refresh : bool, optional
        When ``True``, bypass cached dated-release listings.
    json_output : bool, optional
        When ``True``, output JSON.

    Returns
    -------
    str
        Filtered release catalogue when any selector or ``json_output`` is
        given, otherwise the table of supported species.

    Raises
    ------
    ValueError
        If the selectors conflict.
    """
    from ..dated_releases import AVAILABLE_EXPIRE_AFTER, fetch_all_dated_releases

    if selects_catalogue(
        species=species,
        assembly=assembly,
        assembly_accession=assembly_accession,
        provider=provider,
        reference_name=reference_name,
        include_alt=include_alt,
        json_output=json_output,
    ):
        return available_catalogue(
            species=species,
            assembly=assembly,
            assembly_accession=assembly_accession,
            provider=provider,
            reference_name=reference_name,
            include_alt=include_alt,
            refresh=refresh,
            json_output=json_output,
        )
    return format_available_species(
        dated_releases=fetch_all_dated_releases(
            expire_after=0 if refresh else AVAILABLE_EXPIRE_AFTER,
        )
    )
