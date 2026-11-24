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

"""Shared text formatting, path display, and filesystem summary helpers for the CLI."""

import logging
import os
from pathlib import Path
import sys

logger = logging.getLogger(__name__)

DIVISION_LABELS = (
    ("vertebrates", "Vertebrates"),
    ("metazoa", "Invertebrates"),
    ("plants", "Plants"),
    ("fungi", "Fungi"),
    ("protists", "Protists"),
    ("bacteria", "Bacteria"),
)


def release_sort_key(release):
    """Return a sort key ordering numbered releases before dated releases.

    Parameters
    ----------
    release : int or str
        Release identifier such as ``75``, ``"110"``, or ``"2026_04"``.

    Returns
    -------
    tuple[int, int, str]
        Sort key placing numbered releases in ascending numeric order followed
        by dated or custom release labels in lexicographic order.
    """
    release = str(release)
    return (0, int(release), "") if release.isdigit() else (1, 0, release)


def display_path(path):
    """Shorten paths under the user's home directory with a leading ``~``.

    Parameters
    ----------
    path : str or os.PathLike
        Filesystem path to format.

    Returns
    -------
    str
        Path string with ``~`` substituted for the home prefix when applicable.
    """
    home = os.path.expanduser("~")
    path = os.fspath(path)
    return "~" + path[len(home):] if path == home or path.startswith(home + os.sep) else path


def format_table(header, rows, use_color=None):
    """Render an aligned plain-text table with an optional bold header.

    Cells are converted via ``str``, so integer releases render cleanly.

    Parameters
    ----------
    header : sequence of str
        Column titles.
    rows : sequence of sequence of object
        Data rows matching the column count of ``header``.
    use_color : bool, optional
        Whether to style the header with ANSI bold codes. When ``None``
        (default), enabled only when ``sys.stdout`` is a TTY.

    Returns
    -------
    str
        Formatted multi-line table string.
    """
    if use_color is None:
        use_color = sys.stdout.isatty()
    rows = [[str(cell) for cell in row] for row in rows]
    widths = [max(len(cell) for cell in column) for column in zip(header, *rows)]

    def format_line(cells):
        return "  ".join(cell.ljust(width) for cell, width in zip(cells, widths)).rstrip()

    bold, reset = ("\x1b[1m", "\x1b[0m") if use_color else ("", "")
    return "\n".join([bold + format_line(header) + reset] + [format_line(row) for row in rows])


def format_release_range(start, end):
    """Format an inclusive release range, collapsing single-release ranges.

    Parameters
    ----------
    start : int
        First supported release number.
    end : int
        Last supported release number.

    Returns
    -------
    str
        ``"54"`` when ``start == end``, otherwise ``"<start>–<end>"`` using an
        en-dash.
    """
    if start == end:
        return str(start)
    return "%d–%d" % (start, end)


def species_display_name(species):
    """Return the preferred human-readable common name for ``species``.

    Parameters
    ----------
    species : Species
        Registered species object.

    Returns
    -------
    str
        First registered synonym if available, otherwise ``species.latin_name``.
    """
    if species.synonyms:
        return species.synonyms[0]
    return species.latin_name


def genome_description(genome):
    """Return a concise human-readable description of a genome or release.

    Parameters
    ----------
    genome : Genome
        ``EnsemblRelease``, ``PinnedGenome``, or custom ``Genome`` instance.

    Returns
    -------
    str
        Label such as ``"human GRCh38 release 110"`` or ``"GRCh38 custom 1"``.
    """
    from ..ensembl_release import EnsemblRelease
    from ..pinned_genomes import PinnedGenome

    if isinstance(genome, EnsemblRelease):
        species = genome.species
        name = (
            species.synonyms[0] if species.synonyms else species.latin_name
        ) if species is not None else "assembly"
        return "%s %s release %s" % (name, genome.reference_name, genome.release)
    if isinstance(genome, PinnedGenome):
        species = genome.species
        name = (
            species.synonyms[0] if species.synonyms else species.latin_name
        ) if species is not None else "assembly"
        return "%s %s %s %s" % (
            name,
            genome.reference_name,
            genome.annotation_name,
            genome.annotation_version,
        )
    description = "%s %s" % (genome.reference_name, genome.annotation_name)
    if genome.annotation_version is not None:
        description += " %s" % genome.annotation_version
    return description


def directory_size(path):
    """Sum file and symlink sizes under ``path`` without following links.

    Parameters
    ----------
    path : str or os.PathLike
        Directory to measure.

    Returns
    -------
    int
        Total size in bytes.
    """
    size = 0
    for root, directories, files in os.walk(path):
        for name in files + directories:
            entry = os.path.join(root, name)
            if os.path.isfile(entry) or os.path.islink(entry):
                size += os.lstat(entry).st_size
    return size


def directory_stats(path):
    """Return ``(exists, total_bytes, file_count)`` for ``path`` read-only.

    Does not follow symlinks or create any directories.

    Parameters
    ----------
    path : str or os.PathLike
        Directory path to inspect.

    Returns
    -------
    tuple[bool, int, int]
        ``(exists, total_bytes, file_count)``.
    """
    if not os.path.exists(path):
        return False, 0, 0
    total_bytes = 0
    file_count = 0
    for root, _dirs, files in os.walk(path, followlinks=False):
        for name in files:
            try:
                st = os.lstat(os.path.join(root, name))
            except OSError:
                continue
            total_bytes += st.st_size
            file_count += 1
    return True, total_bytes, file_count


def subdirectories(path):
    """Return sorted non-hidden child directories of ``path``, ignoring errors.

    Parameters
    ----------
    path : str or os.PathLike
        Parent directory to scan.

    Returns
    -------
    list of pathlib.Path
        Sorted visible subdirectories; empty if ``path`` is missing or unreadable.
    """
    try:
        return sorted(
            entry for entry in Path(path).iterdir()
            if not entry.name.startswith(".") and entry.is_dir()
        )
    except FileNotFoundError:
        return []
    except OSError as error:
        logger.warning("Cannot read %s: %s", path, error)
        return []
