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

"""Command-line entry point and action dispatch for ``pyensembl``."""

from contextlib import contextmanager
import os

import argcomplete

from . import (
    available,
    delete,
    inspect as inspect_action,
    install,
    listing,
    logging_config,
    parser as parser_module,
    selection,
    verify,
)


@contextmanager
def offline_override(offline):
    """Temporarily set ``DATACACHE_OFFLINE=1`` for the duration of a CLI call.

    Parameters
    ----------
    offline : bool
        When ``True``, sets ``DATACACHE_OFFLINE=1`` and restores the previous
        value on exit. When ``False``, acts as a no-op.
    """
    if not offline:
        yield
        return
    previous = os.environ.get("DATACACHE_OFFLINE")
    os.environ["DATACACHE_OFFLINE"] = "1"
    try:
        yield
    finally:
        if previous is None:
            os.environ.pop("DATACACHE_OFFLINE", None)
        else:
            os.environ["DATACACHE_OFFLINE"] = previous


def run_action(args):
    """Execute the parsed CLI action described by ``args``.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments from :data:`pyensembl.cli.parser.parser`.
    """
    active_parser = parser_module.parser
    if args.action == "prune":
        try:
            print(delete.run_prune(args.dry_run))
        except ValueError as error:
            active_parser.error(str(error))
        return
    if args.dry_run:
        active_parser.error("--dry-run requires prune")
    if (
        args.action in parser_module.DELETE_ACTIONS
        and not args.release
        and not (
            args.gtf
            or args.transcript_fasta
            or args.protein_fasta
            or (args.genome_fasta_path and args.annotation_name)
        )
    ):
        active_parser.error("%s requires an explicit --release" % args.action)
    if args.action == "dirs":
        print(listing.format_cache_dirs(json_output=args.json))
        return
    if args.action == "list":
        if args.json:
            print(listing.format_installed_genomes_json(check_genome_fasta=args.check_genome_fasta))
        else:
            print(listing.format_installed_genomes(check_genome_fasta=args.check_genome_fasta))
        return
    if args.action == "available":
        from datacache import OfflineError

        available_kwargs = dict(
            species=args.species,
            assembly=args.assembly,
            assembly_accession=args.assembly_accession,
            provider=args.provider,
            reference_name=args.reference_name,
            include_alt=args.include_alt,
            json_output=args.json,
        )
        try:
            print(available.run_available(**available_kwargs, refresh=args.refresh))
        except OfflineError as error:
            active_parser.exit(1, "%s: error: %s\n" % (active_parser.prog, error))
        except ValueError as error:
            # Conflicting selectors are argument errors; failures while
            # reading release listings keep their tracebacks.
            if not available.selects_catalogue(**available_kwargs):
                raise
            active_parser.error(str(error))
        return

    try:
        genomes = selection.collect_selected_genomes(args)
    except ValueError as error:
        active_parser.error(str(error))
    if len(genomes) == 0:
        active_parser.error("No genomes selected")

    if args.action == "verify":
        from datacache import OfflineError

        try:
            verify.require_ensembl_releases(genomes)
        except ValueError as error:
            active_parser.error(str(error))
        # Checksum and manifest errors keep their tracebacks.
        try:
            text, ok = verify.run_verify(genomes, args.check_remote, args.json)
        except OfflineError as error:
            active_parser.exit(1, "%s: error: %s\n" % (active_parser.prog, error))
        print(text)
        if not ok:
            raise SystemExit(1)
    elif args.action == "inspect":
        print(inspect_action.inspect_genomes(genomes, args.check_genome_fasta, args.json))
    elif args.action in parser_module.DELETE_ACTIONS:
        for genome in genomes:
            delete.delete_genome_files(genome, args.action)
    elif args.action == "install":
        from datacache import OfflineError
        from ..dated_releases import UnpublishedDateError

        for genome in genomes:
            try:
                install.run_install(genome, args.only_genome_fasta, args.overwrite)
            except UnpublishedDateError as error:
                active_parser.error(str(error))
            except (install.InstallError, OfflineError) as error:
                active_parser.exit(1, "%s: error: %s\n" % (active_parser.prog, error))
    else:
        raise ValueError("Invalid action: %s" % args.action)


def main(argv=None):
    """Command-line entry point for the ``pyensembl`` console script.

    Parameters
    ----------
    argv : list of str, optional
        Arguments after the program name; defaults to ``sys.argv[1:]``.
    """
    active_parser = parser_module.parser
    # When the shell asks for completions, this prints them and exits.
    argcomplete.autocomplete(active_parser)
    args = active_parser.parse_args(argv)
    logging_config.configure_logging(verbose=args.verbose)
    try:
        parser_module.validate_arguments(
            args.action, args.json, args.refresh, args.check_remote, args.release
        )
    except ValueError as error:
        active_parser.error(str(error))
    with offline_override(args.offline):
        run_action(args)


run = main
