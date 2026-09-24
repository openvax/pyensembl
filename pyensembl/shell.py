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

"""
Manipulate pyensembl's local cache.

    %(prog)s {install, delete-all-files, delete-index-files, list, available} [--release XXX --species human...]

To install particular Ensembl human release(s):
    %(prog)s install --release 75 77

To install particular Ensembl mouse release(s):
    %(prog)s install --release 75 77 --species mouse

To install the newest supported Ensembl release for a reference assembly:
    %(prog)s install --reference-name GRCh37

To delete all downloaded and cached data for a particular Ensembl release:
    %(prog)s delete-all-files --release 75 --species human

To delete everything except the original GTF and FASTA files:
    %(prog)s delete-index-files --release 75

To list all installed genomes:
    %(prog)s list

To list supported species and their Ensembl release ranges:
    %(prog)s available

To install a genome from source files:
    %(prog)s install \
 --reference-name "GRCh38" \
 --gtf URL_OR_PATH \
 --transcript-fasta URL_OR_PATH \
 --protein-fasta URL_OR_PATH
"""

import argparse
import logging.config
from importlib import resources
import os

from .ensembl_release import EnsemblRelease
from .ensembl_versions import MAX_ENSEMBL_RELEASE
from .genome import Genome
from .reference_name import find_species_by_reference, normalize_reference_name
from .species import Species, find_species_by_name
from .version import __version__

logger = logging.getLogger(__name__)


def configure_logging():
    """Apply pyensembl's console logging configuration.

    This is only invoked from the command-line entrypoint (``run``) so that
    merely importing this module never reconfigures the root logger or
    disables any loggers the host application has already created.
    """
    logging.config.fileConfig(
        str(resources.files("pyensembl") / "logging.conf"),
        disable_existing_loggers=False,
    )


parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument(
    "--version", 
    action="version",
    version='%(prog)s {version}'.format(version=__version__)
)

parser.add_argument(
    "--overwrite",
    default=False,
    action="store_true",
    help="Force download and indexing even if files already exist locally",
)


release_group = parser.add_argument_group("Ensembl release options")
release_group.add_argument(
    "--release",
    type=int,
    nargs="+",
    default=[],
    help=(
        "Ensembl release version(s); required for deletion "
        "(install default=newest supported release for --reference-name, "
        "otherwise %d)" % MAX_ENSEMBL_RELEASE
    ),
)

release_group.add_argument(
    "--species",
    default=[],
    nargs="+",
    help=(
        "Which species to download Ensembl data for "
        "(default=inferred from --reference-name, otherwise human)"
    ),
)

release_group.add_argument(
    "--custom-mirror",
    default=None,
    help="URL and directory to use instead of the default Ensembl FTP server",
)

parser.add_argument(
    "--reference-name",
    type=str,
    default=None,
    help=(
        "Reference assembly, e.g. GRCh37 (case-insensitive for Ensembl). "
        "Selects its newest supported release unless --release is given; "
        "with custom source files, names the custom reference."
    ),
)

path_group = parser.add_argument_group("Custom genome options")

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

path_group.add_argument(
    "--shared-prefix",
    default="",
    help="Add this prefix to URLs or paths specified by --gtf, --transcript-fasta, --protein-fasta",
)

parser.add_argument(
    "action",
    type=lambda arg: arg.lower().strip(),
    choices=(
        "install",
        "delete-all-files",
        "delete-index-files",
        "list",
        "available",
    ),
    help=(
        '"install" will download and index any data that is  not '
        'currently downloaded or indexed. "delete-all-files" will delete all data '
        'associated with a genome annotation. "delete-index-files" deletes '
        "all files other than the original GTF and FASTA files for a genome. "
        '"list" will show you all installed Ensembl genomes. '
        '"available" prints every species and the Ensembl release ranges '
        "supported by pyensembl."
    ),
)


def collect_all_installed_ensembl_releases():
    genomes = []
    for species, release in Species.all_species_release_pairs():
        genome = EnsemblRelease(release, species=species)
        if genome.required_local_files_exist():
            genomes.append(genome)
    return sorted(genomes, key=lambda g: (g.species.latin_name, g.release))


def all_combinations_of_ensembl_genomes(args):
    """
    Use all combinations of species and release versions specified by the
    commandline arguments to return a list of EnsemblRelease or Genome objects.
    A reference name constrains the species and releases; omitted values are
    inferred from that reference's newest supported release.
    The results will typically be of type EnsemblRelease unless the
    --custom-mirror argument was given.
    """
    species_list = args.species if args.species else ["human"]
    release_list = args.release if args.release else [MAX_ENSEMBL_RELEASE]
    if args.reference_name is not None:
        reference_name = normalize_reference_name(args.reference_name)
        reference_species = find_species_by_reference(reference_name)
        species_list = args.species or [reference_species.latin_name]
        for species_name in species_list:
            species = find_species_by_name(species_name)
            if species != reference_species:
                raise ValueError(
                    "Reference %s belongs to %s, not %s"
                    % (reference_name, reference_species.latin_name, species.latin_name)
                )
        first_release, last_release = reference_species.reference_assemblies[reference_name]
        release_list = args.release or [last_release]
        for version in release_list:
            if not first_release <= version <= last_release:
                raise ValueError(
                    "Reference %s supports Ensembl releases %d-%d, not --release %d"
                    % (reference_name, first_release, last_release, version)
                )
    genomes = []
    for species in species_list:
        # Otherwise, use Ensembl release information
        for version in release_list:
            ensembl_release = EnsemblRelease(version, species=species)

            if not args.custom_mirror:
                genomes.append(ensembl_release)
            else:
                # if we're using a custom mirror then we expect the provided
                # URL to be a directory with all the same filenames as
                # would be provided by Ensembl
                gtf_url = os.path.join(
                    args.custom_mirror, os.path.basename(ensembl_release.gtf_url)
                )
                transcript_fasta_urls = [
                    os.path.join(
                        args.custom_mirror, os.path.basename(transcript_fasta_url)
                    )
                    for transcript_fasta_url in ensembl_release.transcript_fasta_urls
                ]
                protein_fasta_urls = [
                    os.path.join(
                        args.custom_mirror, os.path.basename(protein_fasta_url)
                    )
                    for protein_fasta_url in ensembl_release.protein_fasta_urls
                ]
                reference_name = ensembl_release.reference_name
                genome = Genome(
                    reference_name=reference_name,
                    annotation_name="ensembl",
                    annotation_version=version,
                    gtf_path_or_url=gtf_url,
                    transcript_fasta_paths_or_urls=transcript_fasta_urls,
                    protein_fasta_paths_or_urls=protein_fasta_urls,
                )
                genomes.append(genome)
    return genomes


def collect_selected_genomes(args):
    # If specific genome source URLs are provided, use those
    if args.gtf or args.transcript_fasta or args.protein_fasta:
        if args.release:
            raise ValueError(
                "An Ensembl release cannot be specified if "
                "specific paths are also given"
            )
        if not args.reference_name:
            raise ValueError("Must specify a reference name")
        if not args.annotation_name:
            raise ValueError("Must specify the name of the annotation source")

        return [
            Genome(
                reference_name=args.reference_name,
                annotation_name=args.annotation_name,
                annotation_version=args.annotation_version,
                gtf_path_or_url=os.path.join(args.shared_prefix, args.gtf),
                transcript_fasta_paths_or_urls=[
                    os.path.join(args.shared_prefix, transcript_fasta)
                    for transcript_fasta in args.transcript_fasta
                ],
                protein_fasta_paths_or_urls=[
                    os.path.join(args.shared_prefix, protein_fasta)
                    for protein_fasta in args.protein_fasta
                ],
            )
        ]
    else:
        return all_combinations_of_ensembl_genomes(args)


_DIVISION_LABELS = (
    ("vertebrates", "Vertebrates"),
    ("metazoa", "Invertebrates"),
    ("plants", "Plants"),
    ("fungi", "Fungi"),
    ("protists", "Protists"),
    ("bacteria", "Bacteria"),
)


def _format_release_range(start, end):
    # Inclusive range; collapse a single-value range to just that integer.
    if start == end:
        return str(start)
    return "%d–%d" % (start, end)  # en-dash


def _species_display_name(species):
    if species.synonyms:
        return species.synonyms[0]
    return species.latin_name


def format_available_species(use_color=None):
    """
    Render the table printed by the "available" CLI action: every registered
    species and its supported Ensembl release ranges, grouped by division.

    When ``use_color`` is ``None`` (the default), ANSI styling is applied if
    stdout is a TTY and suppressed otherwise.
    """
    import sys

    if use_color is None:
        use_color = sys.stdout.isatty()
    BOLD = "\x1b[1m" if use_color else ""
    DIM = "\x1b[2m" if use_color else ""
    RESET = "\x1b[0m" if use_color else ""

    species_by_division = {key: [] for key, _ in _DIVISION_LABELS}
    for latin_name in sorted(Species._latin_names_to_species):
        species = Species._latin_names_to_species[latin_name]
        species_by_division.setdefault(species.division, []).append(species)

    all_species = [s for group in species_by_division.values() for s in group]
    if not all_species:
        return ""

    def _w(values, fallback):
        return max((len(v) for v in values), default=fallback)

    name_w = _w([_species_display_name(s) for s in all_species], 8)
    asm_w = _w(
        [asm for s in all_species for asm in s.reference_assemblies], 8
    )
    rng_w = _w(
        [
            _format_release_range(start, end)
            for s in all_species
            for (start, end) in s.reference_assemblies.values()
        ],
        8,
    )
    latin_w = _w([s.latin_name for s in all_species], 8)

    col_name = max(name_w, len("Species")) + 2
    col_asm = max(asm_w, len("Assembly")) + 2
    col_rng = max(rng_w, len("Releases")) + 2
    col_latin = max(latin_w, len("Latin name"))
    total_w = col_name + col_asm + col_rng + col_latin

    lines = []
    header_row = "%-*s%-*s%-*s%s" % (
        col_name, "Species",
        col_asm, "Assembly",
        col_rng, "Releases",
        "Latin name",
    )
    lines.append("%s%s%s" % (BOLD, header_row, RESET))
    lines.append("─" * total_w)

    for division_key, division_label in _DIVISION_LABELS:
        members = species_by_division.get(division_key, [])
        if not members:
            continue
        members.sort(key=_species_display_name)
        lines.append("")
        section = "── %s " % division_label
        section += "─" * max(2, total_w - len(section))
        lines.append("%s%s%s" % (BOLD, section, RESET))
        for species in members:
            for i, (asm, (start, end)) in enumerate(
                species.reference_assemblies.items()
            ):
                if i == 0:
                    name_cell = _species_display_name(species)
                    latin_cell = "%s%s%s" % (DIM, species.latin_name, RESET)
                else:
                    name_cell = ""
                    latin_cell = ""
                lines.append(
                    "%-*s%-*s%-*s%s"
                    % (
                        col_name, name_cell,
                        col_asm, asm,
                        col_rng, _format_release_range(start, end),
                        latin_cell,
                    )
                )
    return "\n".join(lines)


def format_installed_genomes(genomes, use_color=None):
    """
    Render the table printed by the "list" CLI action: one row per cached
    Ensembl genome with its species common name, assembly, release, index
    status, and cache directory, in the style of the "available" table.

    Genomes whose source files were downloaded but never indexed are marked
    "not indexed" so users know the first query will spend time indexing.
    When no genomes are cached, returns a friendly message instead of an
    empty table.

    When ``use_color`` is ``None`` (the default), ANSI styling is applied if
    stdout is a TTY and suppressed otherwise.
    """
    import sys

    if use_color is None:
        use_color = sys.stdout.isatty()
    BOLD = "\x1b[1m" if use_color else ""
    RESET = "\x1b[0m" if use_color else ""

    if not genomes:
        return (
            "No Ensembl genomes are installed yet.\n"
            "\n"
            "Use `pyensembl install --reference-name <assembly>` to download\n"
            "and index a genome, or `pyensembl available` to see the\n"
            "supported species and their Ensembl release ranges."
        )

    rows = []
    for genome in genomes:
        if isinstance(genome, EnsemblRelease):
            species = _species_display_name(genome.species)
            assembly = genome.reference_name
            release = str(genome.release)
        else:
            species = ""
            assembly = genome.reference_name
            release = (
                "" if genome.annotation_version is None
                else str(genome.annotation_version)
            )
        status = "indexed" if genome.index_files_exist() else "not indexed"
        rows.append(
            (
                species,
                assembly,
                release,
                status,
                genome.download_cache.cache_directory_path,
            )
        )

    def _w(values, fallback):
        return max((len(value) for value in values), default=fallback)

    name_w = _w([row[0] for row in rows], len("Species"))
    asm_w = _w([row[1] for row in rows], len("Assembly"))
    rel_w = _w([row[2] for row in rows], len("Release"))
    status_w = _w([row[3] for row in rows], len("Status"))
    path_w = _w([row[4] for row in rows], len("Path"))

    col_name = max(name_w, len("Species")) + 2
    col_asm = max(asm_w, len("Assembly")) + 2
    col_rel = max(rel_w, len("Release")) + 2
    col_status = max(status_w, len("Status")) + 2
    total_w = col_name + col_asm + col_rel + col_status + path_w

    lines = []
    header_row = "%-*s%-*s%-*s%-*s%s" % (
        col_name, "Species",
        col_asm, "Assembly",
        col_rel, "Release",
        col_status, "Status",
        "Path",
    )
    lines.append("%s%s%s" % (BOLD, header_row, RESET))
    lines.append("─" * total_w)
    for species, assembly, release, status, path in rows:
        lines.append(
            "%-*s%-*s%-*s%-*s%s"
            % (
                col_name, species,
                col_asm, assembly,
                col_rel, release,
                col_status, status,
                path,
            )
        )
    return "\n".join(lines)


def _genome_description(genome):
    if isinstance(genome, EnsemblRelease):
        species = genome.species
        name = species.synonyms[0] if species.synonyms else species.latin_name
        return "%s %s release %d" % (name, genome.reference_name, genome.release)
    description = "%s %s" % (genome.reference_name, genome.annotation_name)
    if genome.annotation_version is not None:
        description += " %s" % genome.annotation_version
    return description


def _directory_size(path):
    """Sum file sizes without following links outside the cache directory."""
    size = 0
    for root, directories, files in os.walk(path):
        for name in files + directories:
            entry = os.path.join(root, name)
            if os.path.isfile(entry) or os.path.islink(entry):
                size += os.lstat(entry).st_size
    return size


def _delete_genome_files(genome, action):
    if action == "delete-index-files":
        deleted = genome.delete_index_files()
    else:
        directory = genome.download_cache.cache_directory_path
        deleted = []
        if os.path.isdir(directory):
            genome.close()
            size = _directory_size(directory)
            genome.download_cache.delete_cache_directory()
            deleted.append((directory, size))
    for path, size in deleted:
        print("Deleted %s (%s bytes)" % (path, format(size, ",")))
    if not deleted:
        print("Nothing to delete for %s" % _genome_description(genome))


def run():
    configure_logging()
    args = parser.parse_args()
    if (
        args.action in ("delete-all-files", "delete-index-files")
        and not args.release
        and not (args.gtf or args.transcript_fasta or args.protein_fasta)
    ):
        parser.error("%s requires an explicit --release" % args.action)
    if args.action == "list":
        # TODO: how do we also identify which non-Ensembl genomes are
        # installed?
        genomes = collect_all_installed_ensembl_releases()
        print(format_installed_genomes(genomes))
    elif args.action == "available":
        print(format_available_species())
    else:
        try:
            genomes = collect_selected_genomes(args)
        except ValueError as error:
            parser.error(str(error))

        if len(genomes) == 0:
            logger.error("ERROR: No genomes selected!")
            parser.print_help()

        for genome in genomes:
            logger.info("Running '%s' for %s", args.action, genome)
            if args.action in ("delete-all-files", "delete-index-files"):
                _delete_genome_files(genome, args.action)
            elif args.action == "install":
                genome.download(overwrite=args.overwrite)
                genome.index(overwrite=args.overwrite)
            else:
                raise ValueError("Invalid action: %s" % args.action)
