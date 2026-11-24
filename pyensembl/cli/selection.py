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

"""Resolve CLI arguments into concrete ``Genome`` / ``EnsemblRelease`` selections."""

from ..ensembl_versions import is_dated_release


def select_pinned_genome(args):
    """Resolve a pinned assembly selection from parsed CLI arguments, or return ``None``.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    list of PinnedGenome or None
        Single-element list containing the selected :class:`PinnedGenome`, or
        ``None`` when the arguments do not name a pinned assembly.

    Raises
    ------
    ValueError
        If assembly selectors, species, provider, release, or mirror flags
        conflict with the pinned assembly definition.
    """
    from ..pinned_genomes import PinnedGenome, lookup_pinned_alias, resolve_pinned_assembly
    from ..species import find_species_by_name

    reference_name = args.reference_name
    assembly = getattr(args, "assembly", None)
    accession = getattr(args, "assembly_accession", None)
    selectors = [value for value in (reference_name, assembly, accession) if value is not None]
    if not any(lookup_pinned_alias(value) is not None for value in selectors):
        return None

    resolved = None
    for value in selectors:
        candidate = resolve_pinned_assembly(value)
        if resolved is not None and candidate.name != resolved.name:
            raise ValueError(
                "Conflicting assembly selectors: %s vs %s"
                % (resolved.reference_name, candidate.reference_name)
            )
        resolved = candidate

    if args.species:
        expected = find_species_by_name(resolved.species)
        for species_name in args.species:
            species = find_species_by_name(species_name)
            if species != expected:
                raise ValueError(
                    "Reference %s belongs to %s, not %s"
                    % (resolved.reference_name, expected.latin_name, species.latin_name)
                )

    provider = getattr(args, "provider", None)
    if provider is not None and provider != resolved.annotation_name:
        raise ValueError(
            "Pinned assembly %s uses provider %s, not %s"
            % (resolved.reference_name, resolved.annotation_name, provider)
        )
    if getattr(args, "include_alt", None) is not None:
        raise ValueError("--include-alt applies only to dated Ensembl Platform releases")
    if getattr(args, "server", None) or args.custom_mirror:
        raise ValueError(
            "Pinned assemblies use checksum-verified URLs and cannot be combined with "
            "--server or --custom-mirror"
        )
    if (args.with_genome_fasta or args.only_genome_fasta) and not args.genome_fasta_path:
        raise ValueError(
            "Pinned assembly %s requires --genome-fasta-path for reference DNA"
            % resolved.reference_name
        )

    for version in args.release:
        if version != "latest" and str(version) != str(resolved.annotation_version):
            raise ValueError(
                "%s has one pinned annotation (%s %s); release %s is not available"
                % (
                    resolved.reference_name,
                    resolved.annotation_name,
                    resolved.annotation_version,
                    version,
                )
            )

    return [
        PinnedGenome(
            resolved.name,
            genome_fasta_path_or_url=args.genome_fasta_path,
            cache_directory_path=getattr(args, "cache_directory", None),
        )
    ]


def all_combinations_of_ensembl_genomes(args):
    """Build all ``EnsemblRelease`` or mirror ``Genome`` objects requested by ``args``.

    A reference name constrains the species and releases; omitted values are
    inferred from that reference's newest supported release. The results are
    typically ``EnsemblRelease`` instances unless ``--custom-mirror`` is given
    or a pinned assembly is selected.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    list of Genome
        Selected ``EnsemblRelease``, ``PinnedGenome``, or mirror ``Genome``
        instances.
    """
    from ..ensembl_release import EnsemblRelease
    from ..reference_name import find_species_by_reference, normalize_reference_name
    from ..species import find_species_by_name

    pinned = select_pinned_genome(args)
    if pinned is not None:
        return pinned

    assembly = getattr(args, "assembly", None)
    accession = getattr(args, "assembly_accession", None)
    if getattr(args, "server", None) and args.custom_mirror:
        raise ValueError("Choose --server for a structured mirror or --custom-mirror for a flat mirror")
    if (assembly or accession) and len(args.species) > 1:
        raise ValueError("An explicit assembly requires at most one --species")
    species_list = args.species or ([None] if assembly or accession else ["human"])
    release_list = args.release
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
            if version == "latest" or is_dated_release(version):
                # The shared resolver handles historical dated assemblies too.
                EnsemblRelease(version, species=reference_species, assembly=reference_name)
            elif not first_release <= version <= last_release:
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
        if not release_list and (assembly or accession):
            raise ValueError("--assembly requires an explicit --release")
        for version in release_list or [find_species_by_name(species).newest_release]:
            selected_assembly = assembly
            if args.reference_name is not None:
                if assembly is not None:
                    probe = EnsemblRelease(
                        version,
                        species=species,
                        assembly=assembly,
                        assembly_accession=accession,
                    )
                    if probe.assembly != reference_name and probe.reference_name != reference_name:
                        raise ValueError("--assembly conflicts with --reference-name")
                else:
                    selected_assembly = reference_name
            if getattr(args, "cache_directory", None) and (
                len(species_list) != 1 or len(release_list or [version]) != 1
            ):
                raise ValueError("--cache-directory requires exactly one dataset")
            ensembl_release = EnsemblRelease(
                version,
                species=species,
                genome_fasta=genome_fasta,
                genome_fasta_type=args.genome_fasta_type,
                genome_fasta_mask=args.masked,
                assembly=selected_assembly,
                assembly_accession=accession,
                provider=getattr(args, "provider", None),
                include_alt=getattr(args, "include_alt", None),
                reference_name=getattr(args, "cache_reference_name", None),
                cache_directory_path=getattr(args, "cache_directory", None),
                server=getattr(args, "server", None),
            )

            if not args.custom_mirror:
                genomes.append(ensembl_release)
            else:
                genomes.append(
                    ensembl_release.with_mirror(
                        args.custom_mirror,
                        genome_fasta_path_or_url=args.genome_fasta_path,
                    )
                )
    return genomes


def collect_selected_genomes(args):
    """Select custom or Ensembl genomes described by parsed CLI arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    list of Genome
        Concrete genomes selected by ``args``.

    Raises
    ------
    ValueError
        If any source path is explicitly empty or options conflict.
    """
    from ..genome import Genome

    # An explicitly empty source is an argument error, not an omitted option.
    # In particular, do not let --gtf "" or --custom-mirror "" select the
    # default human download through the truthiness checks below.
    sources = [
        ("--gtf", args.gtf),
        ("--genome-fasta-path", args.genome_fasta_path),
        ("--custom-mirror", args.custom_mirror),
        ("--server", getattr(args, "server", None)),
        *(("--transcript-fasta", path) for path in args.transcript_fasta),
        *(("--protein-fasta", path) for path in args.protein_fasta),
    ]
    for option, value in sources:
        if value == "":
            raise ValueError("%s must not be empty" % option)

    if args.gtf or args.transcript_fasta or args.protein_fasta or (
        args.genome_fasta_path and args.annotation_name
    ):
        if any(
            getattr(args, key, None) is not None
            for key in (
                "assembly",
                "assembly_accession",
                "provider",
                "include_alt",
                "cache_reference_name",
                "cache_directory",
                "server",
            )
        ):
            raise ValueError("Ensembl assembly options cannot be combined with custom source files")
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
                gtf_path_or_url=args.gtf or None,
                transcript_fasta_paths_or_urls=list(args.transcript_fasta),
                protein_fasta_paths_or_urls=list(args.protein_fasta),
                genome_fasta_path_or_url=args.genome_fasta_path or None,
            )
        ]
    return all_combinations_of_ensembl_genomes(args)
