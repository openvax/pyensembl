"""One offline view of numbered and dated annotations, with scoped evidence."""

from copy import deepcopy
from functools import lru_cache
from importlib.resources import files
import json
from pathlib import Path
import re

from .assembly import (
    _matches, catalogue, check_accession, normalize_species, numbered_assembly,
    species_matches, validate_assembly,
)
from .dated_releases import _cached_dates, available_annotation_dates
from .download_cache import cache_root


@lru_cache(maxsize=1)
def comparison_evidence():
    return json.loads(files("pyensembl").joinpath("release_comparisons.json").read_text())


def _newest_number(species):
    return max(last for first, last in species.reference_assemblies.values())


def available_releases(species=None, *, assembly=None, assembly_accession=None,
                       provider=None, include_alt=None, refresh=False):
    """List metadata for the selected assembly's numbered and dated releases.

    Uses bundled metadata and existing directory listings without network or
    cache writes unless ``refresh=True``. With no assembly selection, use the
    species' latest numbered assembly. Entries are grouped by kind for display;
    their list order does not imply a cross-platform chronology. ``latest``
    marks the conservative, verified choice used by EnsemblRelease("latest").
    Explicit provider/coverage options select the dated platform only.
    """
    species = normalize_species(species, assembly, assembly_accession)
    if assembly_accession is None and isinstance(assembly, str) and re.fullmatch(
        r"GC[AF]_\d{9}\.[1-9]\d*", assembly.strip().upper(),
    ):
        assembly_accession = assembly.strip().upper()
        assembly = assembly_accession
    if assembly_accession is not None:
        check_accession(assembly_accession)
    if provider is not None and (not isinstance(provider, str) or not re.fullmatch(r"[a-z][a-z0-9_]*", provider)):
        raise ValueError("provider must be one lowercase directory name")
    if include_alt is not None and not isinstance(include_alt, bool):
        raise TypeError("include_alt must be bool")
    selector = assembly_accession or assembly
    if selector is None and species is not None:
        selector = species.which_reference(_newest_number(species))
    rows = []
    if species is not None and provider is None and include_alt is None:
        for first, last in species.reference_assemblies.values():
            for release in range(first, last + 1):
                metadata = numbered_assembly(species, release)
                if _matches(metadata, selector):
                    validate_assembly(metadata, assembly, assembly_accession)
                    rows.append(dict(metadata, release=release, kind="numbered",
                                     provider="ensembl", coverage="legacy", subset=False,
                                     comparisons=[], latest=False))
    datasets = [r for r in catalogue()["dated"] if _matches(r, selector)
                and (species is None or species_matches(species, r))
                and (provider is None or provider == r["provider"])]
    if not datasets and assembly_accession:
        # Arbitrary assemblies still use the existing listing/discovery route.
        from .assembly import dated_assembly
        record, _ = dated_assembly("9999_12", species, assembly, assembly_accession, provider)
        datasets = [record]
    for record in datasets:
        validate_assembly(record, assembly, assembly_accession)
        accession, selected_provider = record["assembly_accession"], record["provider"]
        dates = (available_annotation_dates(accession, selected_provider, refresh=True)
                 if refresh else _cached_dates(accession, selected_provider))
        dates = sorted(set(record.get("dates", [])) | set(dates or []))
        # Match ordinary construction's default coverage, including accession
        # selection's primary-only default. No cache directory is constructed.
        alt = (include_alt if include_alt is not None else
               not assembly_accession and not str(assembly or "").upper().startswith(("GCA_", "GCF_"))
               and species is not None and species.dated_reference_name == "GRCh38"
               and tuple(species.dated_releases or ()) == (accession, selected_provider))
        for date in dates:
            comparisons = [r for r in comparison_evidence()["comparisons"]
                           if r["assembly_accession"] == accession
                           and r["provider"] == selected_provider and r["dated"] == date
                           and r["include_alt"] == alt]
            subset = accession == "GCA_000001405.29" and selected_provider == "ensembl" and date == "2024_11"
            rows.append(dict(
                (k, deepcopy(v)) for k, v in record.items() if k != "dates"
            ) | dict(release=date, kind="dated", provider=selected_provider,
                     coverage="including-alt" if alt else "primary", subset=subset,
                     comparisons=deepcopy(comparisons), latest=False))
    if not rows:
        raise ValueError("No known releases for this species/assembly/provider selection; use explicit API discovery or refresh")
    numbered = [r for r in rows if r["kind"] == "numbered"]
    if numbered:
        selected = max(numbered, key=lambda r: r["release"])
        reason = "Newest frozen numbered release; dated datasets do not replace it without verified ordering"
    else:
        eligible = [r for r in rows if not r["subset"]]
        scopes = {(r["assembly_accession"], r["provider"]) for r in eligible}
        if len(scopes) == 1 and eligible:
            selected = max(eligible, key=lambda r: r["release"])
            reason = "Newest known annotation date within one assembly and provider"
        else:
            selected, reason = None, "Choose an assembly accession/provider with a non-subset dataset"
    # A reviewed future ordering can advance the opt-in selector without ever
    # treating a newly discovered date as numerically newer than release 116.
    for rule in comparison_evidence().get("preferred", []):
        matches = [r for r in rows if all(r.get(k) == rule[k] for k in
                   ("release", "assembly_accession", "provider", "coverage")) and not r["subset"]]
        if matches:
            selected, reason = matches[0], rule["reason"]
    for row in rows:
        row["latest"] = row is selected
        row["selection_reason"] = reason if row is selected else None
    return rows


def resolve_latest(**selection):
    rows = available_releases(**selection)
    selected = [row for row in rows if row["latest"]]
    if len(selected) != 1:
        raise ValueError("Latest is ambiguous; specify an assembly accession and provider")
    return selected[0]


def cached_dated_releases():
    """Discover existing standard and explicit-accession dated cache layouts.

    Only directory names are inspected. Construction, source inspection and
    selection remain offline, including when this runs during package import.
    Arbitrary user-selected cache roots still require explicit construction.
    """
    from .ensembl_release import EnsemblRelease
    from .species import Species

    root = Path(cache_root())
    try:
        references = list(root.iterdir())
    except OSError:
        return []
    genomes = []
    for reference in references:
        if not reference.is_dir():
            continue
        try:
            directories = list(reference.iterdir())
        except OSError:
            continue
        for directory in directories:
            match = re.fullmatch(
                r"ensembl(?:-([a-z][a-z0-9_]*)-(GC[AF]_\d{9}\.[1-9]\d*)-(primary|including-alt)-)?"
                r"(\d{4}_(?:0[1-9]|1[0-2]))", directory.name,
            )
            if not match or not directory.is_dir():
                continue
            provider, accession, coverage, date = match.groups()
            if accession:
                options = [dict(assembly_accession=accession, provider=provider,
                                include_alt=coverage == "including-alt",
                                reference_name=reference.name)]
            else:
                options = [dict(species=s) for s in Species._latin_names_to_species.values()
                           if s.dated_releases and s.dated_reference_name == reference.name]
            for selection in options:
                try:
                    genome = EnsemblRelease(date, **selection)
                except ValueError:
                    continue
                if Path(genome.download_cache.cache_directory_path) == directory:
                    genomes.append(genome)
    return genomes
