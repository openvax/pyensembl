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

    %(prog)s {install, delete-all-files, delete-index-files, list, available, prune} [--release XXX --species human...]

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

To list installed genomes, whether they are indexed, and their reference DNA:
    %(prog)s list

To also install reference DNA, or to remove DNA no installed release uses:
    %(prog)s install --release 75 --with-genome-fasta
    %(prog)s prune --dry-run

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
import logging
import os
import sys

from .ensembl_release import EnsemblRelease
from .ensembl_versions import MAX_ENSEMBL_RELEASE
from .genome import Genome
from .genome_fasta import GenomeFasta
from .genome_fasta_cache import (
    canonical_source_options,
    dna_cache_lock,
    dna_cache_root,
    prune_genome_fastas,
)
from .reference_name import find_species_by_reference, normalize_reference_name
from .species import Species, find_species_by_name
from .version import __version__

logger = logging.getLogger(__name__)


class _ConciseFormatter(logging.Formatter):
    """Plain progress messages; warnings and errors keep their level."""

    def format(self, record):
        message = super().format(record)
        if record.levelno >= logging.WARNING:
            return "%s: %s" % (record.levelname.lower(), message)
        return message


_cli_handler = None


def configure_logging(verbose=False):
    """Show progress on stderr: one line per step, or every detail if verbose.

    Only the command-line entrypoint (``run``) calls this, so importing
    pyensembl never changes the host application's logging.
    """
    global _cli_handler
    handler = logging.StreamHandler()  # stderr; stdout carries results
    handler.setFormatter(
        logging.Formatter("%(asctime)s %(name)s %(levelname)s: %(message)s")
        if verbose
        else _ConciseFormatter()
    )
    for name, level in (
        ("pyensembl", logging.DEBUG if verbose else logging.INFO),
        # datacache reports each download step and SQL statement.
        ("datacache", logging.DEBUG if verbose else logging.WARNING),
    ):
        package_logger = logging.getLogger(name)
        if _cli_handler is not None:
            package_logger.removeHandler(_cli_handler)
        package_logger.setLevel(level)
        package_logger.addHandler(handler)
    _cli_handler = handler


parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument(
    "--version", 
    action="version",
    version='%(prog)s {version}'.format(version=__version__)
)

parser.add_argument(
    "-v",
    "--verbose",
    action="store_true",
    help="Show detailed progress, including download and database steps",
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

dna_group = parser.add_argument_group("Optional reference DNA")
dna_group.add_argument("--with-genome-fasta", action="store_true",
                       help="Also download and index reference DNA (several GB for human)")
dna_group.add_argument("--only-genome-fasta", action="store_true",
                       help="Download/index only reference DNA, without annotation or transcript data")
dna_group.add_argument("--genome-fasta-path", default=None,
                       help="Attach a local reference FASTA; custom Genome sources also accept a URL")
dna_group.add_argument("--genome-fasta-type", choices=("toplevel", "primary_assembly"),
                       default="toplevel", help="Toplevel includes patch/haplotype contigs (default)")
dna_group.add_argument("--masked", choices=("none", "soft", "hard"), default="none",
                       help="Masking of downloaded reference DNA (default: none)")
dna_group.add_argument("--check-genome-fasta", action="store_true",
                       help="With list, check existing DNA indexes without downloading or rebuilding")
dna_group.add_argument("--orphan-genome-fastas", action="store_true",
                       help="With prune: the default and only target; accepted for compatibility")
dna_group.add_argument("--dry-run", action="store_true",
                       help="With prune, report candidates without deleting files")

parser.add_argument(
    "action",
    type=lambda arg: arg.lower().strip(),
    choices=(
        "install",
        "delete-all-files",
        "delete-index-files",
        "list",
        "available",
        "prune",
    ),
    help=(
        '"install" will download and index any data that is  not '
        'currently downloaded or indexed. "delete-all-files" will delete all data '
        'associated with a genome annotation. "delete-index-files" deletes '
        "all files other than the original GTF and FASTA files for a genome. "
        '"list" shows installed genomes and whether they are indexed. '
        '"available" prints every species and the Ensembl release ranges '
        'supported by pyensembl. "prune" removes shared reference DNA that '
        "no installed release uses."
    ),
)


def genome_fasta_status(cache_directory, check=False):
    """Describe DNA recorded in a genome cache directory, or None if none.

    A broken reference is reported rather than raised, so one bad release
    cannot hide the others.
    """
    try:
        dna = GenomeFasta.installed_source(cache_directory)
    except ValueError as error:
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


def _has_files(directory):
    try:
        return any(os.scandir(directory))
    except OSError:
        return False


def _annotation_is_indexed(genome):
    return all(os.path.exists(path) for path in genome._annotation_index_paths())


def _annotation_status(genome):
    """'indexed', 'not indexed' (downloads or partial indexes), or None."""
    if _annotation_is_indexed(genome):
        return "indexed"
    dna = genome._genome_fasta.expected_path if genome.requires_genome_fasta else None
    paths = genome._annotation_index_paths() + [
        path for path in genome.required_local_files() if path != dna
    ]
    return "not indexed" if any(os.path.exists(path) for path in paths) else None


def _display_path(path):
    home = os.path.expanduser("~")
    path = os.fspath(path)
    return "~" + path[len(home):] if path == home or path.startswith(home + os.sep) else path


def _format_table(header, rows, use_color=None):
    if use_color is None:
        use_color = sys.stdout.isatty()
    widths = [max(len(cell) for cell in column) for column in zip(header, *rows)]

    def line(cells):
        return "  ".join(cell.ljust(width) for cell, width in zip(cells, widths)).rstrip()

    bold, reset = ("\x1b[1m", "\x1b[0m") if use_color else ("", "")
    return "\n".join([bold + line(header) + reset] + [line(row) for row in rows])


def collect_all_installed_ensembl_releases():
    """Ensembl releases with files in the local cache, indexed or not."""
    genomes = [
        EnsemblRelease(release, species=species)
        for species, release in Species.all_species_release_pairs()
    ]
    return sorted(
        (g for g in genomes if _has_files(g.download_cache.cache_directory_path)),
        key=lambda g: (g.species.latin_name, g.release),
    )


def format_installed_genomes(check_genome_fasta=False, use_color=None):
    """A table of genomes in the local cache, or a note that there are none."""
    rows = []
    ensembl_directories = set()
    for genome in collect_all_installed_ensembl_releases():
        directory = genome.download_cache.cache_directory_path
        ensembl_directories.add(os.path.normpath(directory))
        rows.append((
            _species_display_name(genome.species),
            genome.reference_name,
            str(genome.release),
            _annotation_status(genome) or "-",
            genome_fasta_status(directory, check=check_genome_fasta) or "-",
            _display_path(directory),
        ))
    rows.sort(key=lambda row: (row[0], row[1], int(row[2])))
    # Custom genomes (e.g. install --gtf ...) live beside Ensembl releases.
    root = dna_cache_root().parent
    references = sorted(root.iterdir()) if root.is_dir() else []
    for reference in references:
        if reference.name == "dna_cache" or not reference.is_dir():
            continue
        for directory in sorted(reference.iterdir()):
            if (
                os.path.normpath(directory) in ensembl_directories
                or not directory.is_dir()
                or not _has_files(directory)
            ):
                continue
            indexed = any(directory.glob("*.db")) or any(directory.glob("*.pickle"))
            rows.append((
                "custom",
                reference.name,
                directory.name,
                "indexed" if indexed else "not indexed",
                genome_fasta_status(directory, check=check_genome_fasta) or "-",
                _display_path(directory),
            ))
    if not rows:
        return "No genomes installed in %s" % _display_path(root)
    header = ("Species", "Assembly", "Release", "Annotation", "Reference DNA", "Location")
    return _format_table(header, rows, use_color)


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
    dna_requested = args.with_genome_fasta or args.only_genome_fasta
    if args.genome_fasta_path is not None and not args.custom_mirror:
        if not args.genome_fasta_path or "://" in args.genome_fasta_path:
            raise ValueError(
                "--genome-fasta-path must be a local file for an Ensembl release; "
                "use --annotation-name to install a custom genome with a DNA URL"
            )
        genome_fasta = args.genome_fasta_path
    else:
        # A mirror serves Ensembl's filenames; the release only names them.
        genome_fasta = dna_requested
    genomes = []
    for species in species_list:
        # Otherwise, use Ensembl release information
        for version in release_list:
            ensembl_release = EnsemblRelease(
                version, species=species,
                genome_fasta=genome_fasta,
                genome_fasta_type=args.genome_fasta_type,
                genome_fasta_mask=args.masked,
            )

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
                genome_fasta_source = args.genome_fasta_path
                if genome_fasta_source is None and ensembl_release.genome_fasta_urls:
                    genome_fasta_source = os.path.join(
                        args.custom_mirror, os.path.basename(ensembl_release.genome_fasta_urls[0])
                    )
                genome = Genome(
                    reference_name=reference_name,
                    annotation_name="ensembl",
                    annotation_version=version,
                    gtf_path_or_url=gtf_url,
                    transcript_fasta_paths_or_urls=transcript_fasta_urls,
                    protein_fasta_paths_or_urls=protein_fasta_urls,
                    genome_fasta_path_or_url=genome_fasta_source,
                )
                genomes.append(genome)
    return genomes


def collect_selected_genomes(args):
    # If specific genome source URLs are provided, use those
    if args.gtf or args.transcript_fasta or args.protein_fasta or (
        args.genome_fasta_path and args.annotation_name
    ):
        if args.release:
            raise ValueError(
                "An Ensembl release cannot be specified if "
                "specific paths are also given"
            )
        if not args.reference_name:
            raise ValueError("Must specify a reference name")
        if not args.annotation_name:
            raise ValueError("Must specify the name of the annotation source")
        if (args.with_genome_fasta or args.only_genome_fasta) and not args.genome_fasta_path:
            raise ValueError("Custom genomes require --genome-fasta-path for reference DNA")

        return [
            Genome(
                reference_name=args.reference_name,
                annotation_name=args.annotation_name,
                annotation_version=args.annotation_version,
                gtf_path_or_url=os.path.join(args.shared_prefix, args.gtf) if args.gtf else None,
                transcript_fasta_paths_or_urls=[
                    os.path.join(args.shared_prefix, transcript_fasta)
                    for transcript_fasta in args.transcript_fasta
                ],
                protein_fasta_paths_or_urls=[
                    os.path.join(args.shared_prefix, protein_fasta)
                    for protein_fasta in args.protein_fasta
                ],
                genome_fasta_path_or_url=(
                    os.path.join(args.shared_prefix, args.genome_fasta_path)
                    if args.genome_fasta_path else None
                ),
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


def _install(genome, only_genome_fasta=False, overwrite=False):
    description = _genome_description(genome)
    installed = (
        not genome.requires_genome_fasta or genome._genome_fasta.status() == "indexed"
    ) and (only_genome_fasta or _annotation_is_indexed(genome))
    if installed and not overwrite:
        # Still run the idempotent steps below: they finish partial work.
        logger.info("%s is already installed", description)
    else:
        logger.info("Installing %s", description)
    if only_genome_fasta:
        genome.download_genome_fasta(overwrite=overwrite)
        genome.index_genome_fasta(overwrite=overwrite)
    else:
        genome.download(overwrite=overwrite)
        genome.index(overwrite=overwrite)


def _delete_genome_files(genome, action):
    if action == "delete-index-files":
        deleted = genome.delete_index_files()
    else:
        directory = genome.download_cache.cache_directory_path
        deleted = []
        # Serialize reference removal with shared-DNA install/prune. Shared
        # objects themselves are removed only by the explicit prune action.
        if os.path.isdir(directory):
            with dna_cache_lock():
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
    args = parser.parse_args()
    configure_logging(verbose=args.verbose)
    if args.action == "prune":
        try:
            candidates = prune_genome_fastas(dry_run=args.dry_run)
        except ValueError as error:
            parser.error(str(error))
        for path, size in candidates:
            print("%s %s (%s bytes)" % (
                "Would delete" if args.dry_run else "Deleted", path, format(size, ",")
            ))
        if not candidates:
            print("No orphan genome FASTAs")
        return
    if args.orphan_genome_fastas or args.dry_run:
        parser.error("--orphan-genome-fastas and --dry-run require prune")
    if (
        args.action in ("delete-all-files", "delete-index-files")
        and not args.release
        and not (args.gtf or args.transcript_fasta or args.protein_fasta or
                 (args.genome_fasta_path and args.annotation_name))
    ):
        parser.error("%s requires an explicit --release" % args.action)
    if args.action == "list":
        print(format_installed_genomes(check_genome_fasta=args.check_genome_fasta))
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
            if args.action in ("delete-all-files", "delete-index-files"):
                _delete_genome_files(genome, args.action)
            elif args.action == "install":
                _install(genome, args.only_genome_fasta, args.overwrite)
            else:
                raise ValueError("Invalid action: %s" % args.action)
