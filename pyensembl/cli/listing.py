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

"""Implementation of the ``pyensembl list`` and ``pyensembl dirs`` CLI actions."""

import json
import logging
import os
from pathlib import Path
import re

from ..ensembl_versions import is_dated_release
from .formatting import (
    directory_stats,
    display_path,
    format_table,
    release_sort_key,
    species_display_name,
    subdirectories,
)

logger = logging.getLogger(__name__)

# Reference DNA bookkeeping entries in a genome's cache directory.
DNA_ENTRIES = frozenset({"genome_fasta", "genome_fasta.json", "genome_fasta_refs"})


def dna_status(cache_directory, check=False):
    """Describe reference DNA recorded in a genome cache directory, or ``None``.

    A broken reference record is reported rather than raised, so one bad
    release cannot hide the others.

    Parameters
    ----------
    cache_directory : str or os.PathLike
        Genome cache directory to inspect.
    check : bool, optional
        When ``True``, verify that the FASTA index opens and matches sequence
        lengths on disk.

    Returns
    -------
    str or None
        Human-readable DNA status summary, or ``None`` if no DNA is recorded.
    """
    from ..genome_fasta import GenomeFasta
    from ..genome_fasta_cache import canonical_source_options

    try:
        dna = GenomeFasta.installed_source(cache_directory)
    except (OSError, ValueError) as error:
        logger.warning("Invalid reference DNA record in %s: %s", cache_directory, error)
        return "invalid reference"
    if dna is None:
        return None
    if not dna.remote:
        kind = "local " + os.path.basename(dna.source)
    else:
        options = canonical_source_options(dna.source)
        kind = "downloaded"
        if options is not None:
            fasta_type, mask = options
            kind = fasta_type if mask == "none" else "%s, %s-masked" % (fasta_type, mask)
    return "%s, %s" % (kind, dna.status(check=check))


def genome_fasta_status(genome, check=False):
    """Describe the reference DNA recorded for ``genome``'s cache, or ``None``.

    Parameters
    ----------
    genome : Genome
        Genome instance whose cache directory is inspected.
    check : bool, optional
        When ``True``, verify DNA index integrity.

    Returns
    -------
    str or None
        Reference DNA status string, or ``None`` if none is recorded.
    """
    return dna_status(genome.download_cache.cache_directory_path, check=check)


def directory_annotation_status(directory):
    """Return best-effort annotation status for an unrecognized cache directory.

    Parameters
    ----------
    directory : str or os.PathLike
        Genome cache subdirectory to inspect.

    Returns
    -------
    str or None
        ``"indexed"``, ``"not indexed"``, or ``None`` if no annotation files
        are present.
    """
    from ..database import is_complete_database
    from ..fasta_cache import is_fasta_index_filename

    try:
        names = [
            name for name in os.listdir(directory)
            if not name.startswith(".") and name not in DNA_ENTRIES
        ]
    except OSError:
        return None
    if not names:
        return None
    databases = [name for name in names if name.endswith(".db")]
    if databases:
        complete = all(
            is_complete_database(os.path.join(directory, name)) for name in databases
        )
        return "indexed" if complete else "not indexed"
    if any(name.endswith((".gtf", ".gtf.gz")) for name in names):
        return "not indexed"
    return "indexed" if any(is_fasta_index_filename(name) for name in names) else "not indexed"


def other_genome_labels(reference, annotation):
    """Infer ``(species, release)`` table labels for an unmatched cache directory.

    Parameters
    ----------
    reference : str
        Reference assembly directory name.
    annotation : str
        Annotation subdirectory name (e.g. ``"ensembl2026_04"`` or ``"custom"``).

    Returns
    -------
    tuple[str, str]
        ``(species_label, release_label)`` pair for CLI table display.
    """
    from ..reference_name import find_species_by_reference

    release = annotation[len("ensembl"):] if annotation.startswith("ensembl") else ""
    if not (re.fullmatch(r"[0-9]+", release) or is_dated_release(release)):
        return "custom", annotation
    try:
        species = species_display_name(find_species_by_reference(reference))
    except (KeyError, ValueError):
        species = "unknown"
    return species, release


def cached_dated_releases():
    """Return dated ``EnsemblRelease`` instances that have local cache directories.

    Returns
    -------
    list of EnsemblRelease
        Cached dated release instances discovered on disk.
    """
    from ..release_catalogue import cached_dated_releases as load_cached_dated_releases

    return load_cached_dated_releases()


def ensembl_releases():
    """Return every numbered release plus cached dated and pinned releases.

    Returns
    -------
    list of Genome
        All known numbered ``EnsemblRelease`` objects combined with locally
        cached dated releases and pinned genomes.
    """
    from ..ensembl_release import EnsemblRelease
    from ..pinned_genomes import cached_pinned_genomes
    from ..species import Species

    numbered = [
        EnsemblRelease(release, species=species)
        for species, release in Species.all_species_release_pairs()
    ]
    return numbered + cached_dated_releases() + cached_pinned_genomes()


def has_genome_data(annotation, dna):
    """Return whether a cache entry contains annotation data or reference DNA.

    Parameters
    ----------
    annotation : str or None
        Annotation status string.
    dna : str or None
        Reference DNA status string.

    Returns
    -------
    bool
        ``True`` if either annotation files or reference DNA are present.
    """
    return annotation not in (None, "missing") or dna is not None


def collect_all_installed_ensembl_releases():
    """Return all Ensembl releases with annotation data or reference DNA in cache.

    Returns
    -------
    list of EnsemblRelease
        Installed Ensembl/pinned releases sorted by species latin name and
        release order.
    """
    return sorted(
        (
            genome for genome in ensembl_releases()
            if has_genome_data(genome.annotation_status(), genome_fasta_status(genome))
        ),
        key=lambda genome: (genome.species.latin_name, release_sort_key(genome.release)),
    )


def installed_genome_rows(check_genome_fasta=False):
    """Collect installed genomes as structured row dictionaries with absolute paths.

    Each dictionary contains keys ``species``, ``assembly``, ``release``,
    ``annotation``, ``reference_dna``, and ``cache_directory``.

    Parameters
    ----------
    check_genome_fasta : bool, optional
        When ``True``, verify reference DNA index integrity on disk.

    Returns
    -------
    tuple[list of dict, pathlib.Path]
        ``(rows, cache_root_path)`` sorted deterministically by species,
        assembly, and release.
    """
    from ..download_cache import NON_GENOME_CACHE_DIRECTORIES, cache_root

    rows = []
    ensembl_directories = set()
    ensembl_named = set()
    for genome in ensembl_releases():
        directory = genome.download_cache.cache_directory_path
        ensembl_named.add(os.path.normpath(directory))
        annotation = genome.annotation_status()
        dna = dna_status(directory, check=check_genome_fasta)
        if not has_genome_data(annotation, dna):
            continue  # Possibly a custom genome with the same name; listed below.
        ensembl_directories.add(os.path.normpath(directory))
        rows.append({
            "species": species_display_name(genome.species),
            "assembly": genome.reference_name,
            "release": str(genome.release),
            "annotation": annotation if annotation not in (None, "missing") else None,
            "reference_dna": dna or None,
            "cache_directory": os.fspath(directory),
        })
    # Custom genomes (e.g. install --gtf ...) live beside Ensembl releases.
    root = Path(cache_root())
    for reference in subdirectories(root):
        if reference.name in NON_GENOME_CACHE_DIRECTORIES:
            continue
        for directory in subdirectories(reference):
            if os.path.normpath(directory) in ensembl_directories:
                continue
            annotation = directory_annotation_status(directory)
            dna = dna_status(directory, check=check_genome_fasta)
            if not has_genome_data(annotation, dna):
                continue
            if os.path.normpath(directory) in ensembl_named:
                # Named like an Ensembl release, but holding other files.
                species, release = "custom", directory.name
            else:
                species, release = other_genome_labels(reference.name, directory.name)
            rows.append({
                "species": species,
                "assembly": reference.name,
                "release": release,
                "annotation": annotation or None,
                "reference_dna": dna or None,
                "cache_directory": os.fspath(directory),
            })
    rows.sort(key=lambda r: (r["species"], r["assembly"], release_sort_key(r["release"])))
    return rows, root


def format_installed_genomes(check_genome_fasta=False, use_color=None):
    """Render a human-readable table of genomes in the local cache.

    Parameters
    ----------
    check_genome_fasta : bool, optional
        When ``True``, verify DNA indexes on disk.
    use_color : bool, optional
        Override ANSI header and status styling (defaults to
        :func:`pyensembl.cli.help.color_enabled` on ``sys.stdout``).

    Returns
    -------
    str
        Formatted table string, or a note that no genomes are installed.
    """
    rows, root = installed_genome_rows(check_genome_fasta=check_genome_fasta)
    if not rows:
        return "No genomes installed in %s" % display_path(root)
    header = ("Species", "Assembly", "Release", "Annotation", "Reference DNA", "Location")
    table_rows = [
        (
            r["species"],
            r["assembly"],
            r["release"],
            r["annotation"] or "-",
            r["reference_dna"] or "-",
            display_path(r["cache_directory"]),
        )
        for r in rows
    ]
    return format_table(
        header,
        table_rows,
        use_color,
        status_columns=("Annotation", "Reference DNA"),
    )


def format_installed_genomes_json(check_genome_fasta=False):
    """Serialize the installed genome inventory as formatted JSON.

    Parameters
    ----------
    check_genome_fasta : bool, optional
        When ``True``, verify DNA indexes on disk.

    Returns
    -------
    str
        JSON array string with keys ``species``, ``assembly``, ``release``,
        ``annotation``, ``reference_dna``, and ``cache_directory``.
    """
    rows, _ = installed_genome_rows(check_genome_fasta=check_genome_fasta)
    return json.dumps(rows, indent=2)


def format_cache_dirs(json_output=False):
    """Describe where PyEnsembl's cache lives, why it resolved there, and disk usage.

    Performs read-only inspection only; never creates directories or accesses
    the network.

    Parameters
    ----------
    json_output : bool, optional
        When ``True``, return a JSON document instead of human-readable text.

    Returns
    -------
    str
        Formatted text or JSON summary of the resolved cache root.
    """
    import datacache
    from ..cache_location import resolve_cache_location

    loc = resolve_cache_location()
    exists, total_bytes, file_count = directory_stats(loc.path)
    offline = datacache.offline_source()
    entries = []
    if os.path.isdir(loc.path):
        try:
            children = sorted(os.scandir(loc.path), key=lambda e: e.name)
        except OSError:
            children = []
        for entry in children:
            if not entry.is_dir(follow_symlinks=False):
                continue
            _, sub_bytes, sub_files = directory_stats(entry.path)
            entries.append(dict(
                name=entry.name,
                path=entry.path,
                size_bytes=sub_bytes,
                file_count=sub_files,
            ))
    payload = dict(
        cache_root=loc.path,
        origin=loc.origin,
        detail=loc.detail,
        exists=exists,
        size_bytes=total_bytes,
        file_count=file_count,
        offline=offline is not None,
        offline_source=offline,
        subdirectories=entries,
    )
    if json_output:
        return json.dumps(payload, indent=2)
    lines = [
        "Cache root: %s" % loc.path,
        "Origin:     %s (%s)" % (loc.detail, loc.origin),
        "Size:       %s (%d file(s))" % (datacache.format_bytes(total_bytes), file_count)
        if exists else "Size:       not created yet",
    ]
    if offline:
        lines.append("Offline:    enabled (%s)" % offline)
    if entries:
        lines.append("")
        lines.append(format_table(
            ("Directory", "Size", "Files", "Path"),
            [
                (e["name"], datacache.format_bytes(e["size_bytes"]), str(e["file_count"]), e["path"])
                for e in entries
            ],
        ))
    return "\n".join(lines)
