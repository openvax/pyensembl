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

"""Argument parser, tab-completion helpers, and flag validation for the CLI."""

import argparse
import os
import re

from ..ensembl_versions import normalize_release
from ..version import __version__
from .formatting import release_sort_key
from .help import HelpArgumentParser

DELETE_ACTIONS = ("delete-all-files", "delete-index-files")

ACTIONS_AND_EXAMPLES = """\
actions:
  install             download and index genome data (skips what is already done)
  list                show installed genomes, whether they are indexed, and their DNA
  dirs                show where PyEnsembl's cache lives, why it resolved there, and disk usage
  inspect             inspect selected source files, indexes, and download provenance offline
  verify              verify recorded dated-annotation checksums without repair
  available           show supported species, assemblies, and Ensembl releases
  delete-index-files  delete indexes, keeping downloaded files (needs --release)
  delete-all-files    delete all of a genome's local data (needs --release)
  prune               delete shared reference DNA that no installed release uses

examples:
  pyensembl install --release 75 77                     human releases 75 and 77
  pyensembl install --release 110 --species mouse       a mouse release
  pyensembl install --release 2026_04                   a dated release (new Ensembl platform)
  pyensembl install --reference-name GRCh37             newest release for GRCh37
  pyensembl install --release 110 --with-genome-fasta   also install reference DNA
  pyensembl install --reference-name GRCh38 --annotation-name my_genes \\
      --gtf URL_OR_PATH --transcript-fasta URL_OR_PATH  a custom genome
  pyensembl list
  pyensembl dirs --json
  pyensembl inspect --release 93 --json
  pyensembl delete-all-files --release 75
  pyensembl prune --dry-run
"""


def parse_release_argument(value):
    """Validate and normalize a ``--release`` CLI argument.

    Parameters
    ----------
    value : str
        Raw CLI token (e.g. ``"110"``, ``"2026_04"``, or ``"latest"``).

    Returns
    -------
    int or str
        Normalized release integer, ``YYYY_MM`` date string, or ``"latest"``.

    Raises
    ------
    argparse.ArgumentTypeError
        If ``value`` is not a valid release identifier.
    """
    if value.strip() == "latest":
        return "latest"
    try:
        return normalize_release(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error))


def installed_releases(root):
    """List the releases that have a cache directory, by directory name only.

    Genomes live in ``<root>/<reference>/ensembl<release>``; dated releases
    of non-default datasets use
    ``ensembl-<provider>-<accession>-<coverage>-<date>``. Only directory names
    are read: no file is opened, indexed or downloaded, so this is fast enough
    for tab completion on a large cache.

    Parameters
    ----------
    root : str
        The pyensembl cache directory, e.g. ``cache_root()``.

    Returns
    -------
    list of str
        Distinct release labels such as ``"75"`` or ``"2026_04"``: numbered
        releases in order, then annotation dates in order. Empty if ``root``
        does not exist or cannot be read.
    """
    releases = set()
    try:
        references = [entry for entry in os.scandir(root) if entry.is_dir()]
    except OSError:
        return []
    for reference in references:
        try:
            names = [entry.name for entry in os.scandir(reference.path) if entry.is_dir()]
        except OSError:
            continue
        for name in names:
            match = re.fullmatch(r"ensembl(?:-rapid|-.+-)?(\d+|\d{4}_(?:0[1-9]|1[0-2]))", name)
            if match:
                releases.add(match.group(1))
    return sorted(releases, key=release_sort_key)


def complete_release(prefix, **kwargs):
    """Suggest ``--release`` values: the releases installed in the local cache.

    An argcomplete completer. It lists cache directory names (see
    :func:`installed_releases`) and never touches the network.

    Parameters
    ----------
    prefix : str
        What has been typed so far; argcomplete filters the suggestions by it.
    **kwargs
        Additional keyword arguments passed by argcomplete; unused.

    Returns
    -------
    list of str
        Release labels, e.g. ``["75", "93", "2026_04"]``.
    """
    from ..cache_location import cache_root

    return installed_releases(cache_root())


def complete_species(prefix, **kwargs):
    """Suggest ``--species`` values: every registered species name.

    An argcomplete completer. It reads pyensembl's built-in species registry;
    it never touches the cache or the network.

    Parameters
    ----------
    prefix : str
        What has been typed so far; argcomplete filters the suggestions by it.
    **kwargs
        Additional keyword arguments passed by argcomplete; unused.

    Returns
    -------
    list of str
        Common names (e.g. ``"human"``) and latin names
        (e.g. ``"homo_sapiens"``), sorted.
    """
    from ..species import Species, find_species_by_name

    names = set()
    for latin_name in Species.all_registered_latin_names():
        names.add(latin_name)
        names.update(find_species_by_name(latin_name).synonyms)
    return sorted(names)


def build_parser():
    """Construct and return the ``pyensembl`` CLI :class:`argparse.ArgumentParser`.

    Returns
    -------
    argparse.ArgumentParser
        Configured parser with all actions, groups, and argcomplete completers.
    """
    argument_parser = HelpArgumentParser(
        prog="pyensembl",
        usage="%(prog)s ACTION [options]",
        description="Install and manage the genome data PyEnsembl uses.",
        epilog=ACTIONS_AND_EXAMPLES,
    )

    argument_parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )

    argument_parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Show detailed progress, including download and database steps",
    )

    argument_parser.add_argument(
        "--json",
        action="store_true",
        help="With list, dirs, inspect, available, or verify, print JSON",
    )

    argument_parser.add_argument(
        "--offline",
        action="store_true",
        help="Forbid remote downloads (equivalent to OPENVAX_OFFLINE=1)",
    )

    argument_parser.add_argument(
        "--overwrite",
        default=False,
        action="store_true",
        help="Force download and indexing even if files already exist locally",
    )

    release_group = argument_parser.add_argument_group("Ensembl release options")
    release_group.add_argument(
        "--release",
        type=parse_release_argument,
        nargs="+",
        default=[],
        help=(
            "Ensembl release version(s): numbered, YYYY_MM, or opt-in latest on "
            "the new Ensembl platform; required for deletion "
            "(install default=newest supported release for --reference-name, "
            "otherwise for each species)"
        ),
    ).completer = complete_release

    release_group.add_argument(
        "--species",
        default=[],
        nargs="+",
        help=(
            "Which species to download Ensembl data for "
            "(default=inferred from --reference-name, otherwise human)"
        ),
    ).completer = complete_species

    release_group.add_argument(
        "--assembly",
        help="Optional build, full assembly name, or versioned GCA/GCF accession",
    )
    release_group.add_argument(
        "--assembly-accession",
        help="Explicit accession (uses the compatible accession cache layout)",
    )
    release_group.add_argument("--provider", help="Annotation provider for a dated release")
    release_group.add_argument(
        "--server",
        help="Mirror root preserving Ensembl's directory structure (advanced)",
    )
    release_group.add_argument(
        "--include-alt",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Select dated annotation coverage explicitly",
    )
    release_group.add_argument(
        "--cache-reference-name",
        help="Preserve a custom reference label when reopening an accession dataset",
    )
    release_group.add_argument(
        "--cache-directory",
        help="Explicit existing cache directory for one selected dataset",
    )

    release_group.add_argument(
        "--custom-mirror",
        default=None,
        help=(
            "Flat mirror URL/directory; repeat it to reuse mirror data. In Python use "
            "EnsemblRelease(...).with_mirror(URL). New mirror caches are separate "
            "from default EnsemblRelease caches."
        ),
    )

    argument_parser.add_argument(
        "--reference-name",
        type=str,
        default=None,
        help=(
            "Reference assembly, e.g. GRCh37 (case-insensitive for Ensembl). "
            "Selects its newest supported release unless --release is given; "
            "with custom source files, names the custom reference."
        ),
    )

    path_group = argument_parser.add_argument_group("Custom genome options")

    path_group.add_argument(
        "--annotation-name", default=None, help="Name of annotation source (e.g. refseq)"
    )

    path_group.add_argument(
        "--annotation-version", default=None, help="Version of annotation database"
    )

    path_group.add_argument(
        "--gtf",
        type=str,
        default=None,
        help="URL or local path to a GTF file containing annotations.",
    )

    path_group.add_argument(
        "--transcript-fasta",
        type=str,
        action="append",
        default=[],
        help="URL or local path to a FASTA files containing the transcript "
        "data. This option can be specified multiple times for multiple "
        "FASTA files.",
    )

    path_group.add_argument(
        "--protein-fasta",
        type=str,
        default=[],
        action="append",
        help="URL or local path to a FASTA file containing protein data.",
    )

    dna_group = argument_parser.add_argument_group("Optional reference DNA")
    dna_group.add_argument(
        "--with-genome-fasta",
        action="store_true",
        help="Also download and index reference DNA (several GB for human)",
    )
    dna_group.add_argument(
        "--only-genome-fasta",
        action="store_true",
        help="Download/index only reference DNA, without annotation or transcript data",
    )
    dna_group.add_argument(
        "--genome-fasta-path",
        default=None,
        help="Attach a local reference FASTA; custom Genome sources also accept a URL",
    )
    dna_group.add_argument(
        "--genome-fasta-type",
        choices=("toplevel", "primary_assembly"),
        default="toplevel",
        help="Toplevel includes patch/haplotype contigs (default)",
    )
    dna_group.add_argument(
        "--masked",
        choices=("none", "soft", "hard"),
        default="none",
        help="Masking of downloaded reference DNA (default: none)",
    )
    dna_group.add_argument(
        "--check-genome-fasta",
        action="store_true",
        help="With list or inspect, check DNA indexes without downloading or rebuilding",
    )
    dna_group.add_argument(
        "--dry-run",
        action="store_true",
        help="With prune, report candidates without deleting files",
    )

    argument_parser.add_argument(
        "action",
        type=lambda arg: arg.lower().strip(),
        choices=(
            "install",
            "list",
            "dirs",
            "inspect",
            "verify",
            "available",
            "delete-index-files",
            "delete-all-files",
            "prune",
        ),
        metavar="ACTION",
        help="one of the actions listed below",
    )

    argument_parser.add_argument(
        "--refresh",
        action="store_true",
        help="With available, refresh dated listings from Ensembl",
    )
    argument_parser.add_argument(
        "--check-remote",
        action="store_true",
        help="With verify, compare recorded MD5s with current upstream manifests",
    )
    return argument_parser


parser = build_parser()


def validate_arguments(action, json_output, refresh, check_remote, releases):
    """Check that action-specific flags are used with the right action.

    Parameters
    ----------
    action : str
        The CLI action, e.g. ``"install"``.
    json_output : bool
        Whether ``--json`` was given.
    refresh : bool
        Whether ``--refresh`` was given.
    check_remote : bool
        Whether ``--check-remote`` was given.
    releases : list
        Values of ``--release`` (``"latest"`` or normalized releases).

    Raises
    ------
    ValueError
        With the usage message when a flag does not fit the action.
    """
    if json_output and action not in ("list", "dirs", "inspect", "available", "verify"):
        raise ValueError("--json requires list, dirs, inspect, available, or verify")
    if refresh and action != "available":
        raise ValueError("--refresh requires available")
    if check_remote and action != "verify":
        raise ValueError("--check-remote requires verify")
    if action in DELETE_ACTIONS and "latest" in releases:
        raise ValueError("Deletion requires a concrete --release, not latest")
