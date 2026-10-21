"""Offline assembly identity and selection.

The bundled catalogue describes datasets, not cache directories. In particular,
an upstream assembly name must never replace an installed reference_name.
"""

from functools import lru_cache
from importlib.resources import files
import json
import re

from .species import Species, check_species_object


_ACCESSION = re.compile(r"GC[AF]_\d{9}\.[1-9]\d*\Z")


def check_accession(value):
    if not isinstance(value, str) or not _ACCESSION.fullmatch(value):
        raise ValueError("Expected a versioned GCA/GCF assembly accession")
    return value


@lru_cache(maxsize=1)
def catalogue():
    return json.loads(files("pyensembl").joinpath("assembly_catalog.json").read_text())


def numbered_assembly(species, release):
    record = catalogue()["numbered"].get(species.latin_name, {}).get(str(release))
    if record is not None:
        return dict(record)
    # User-registered species and historical entries without an authoritative
    # accession remain usable offline. Never manufacture an exact name/accession.
    return dict(
        species=species.latin_name, assembly=species.which_reference(release),
        assembly_name=None, assembly_accession=None,
    )


def _matches(record, selector):
    if selector is None:
        return True
    if not isinstance(selector, str) or not selector.strip():
        raise ValueError("assembly must be a build name, full name, or versioned accession")
    value = selector.strip().casefold()
    names = [record.get(k) for k in (
        "assembly", "assembly_name", "assembly_accession", "reference_name",
    )]
    return value in {s.casefold() for s in names + record.get("aliases", []) if s}


def species_matches(species, record):
    """A general species can select a known breed; a different breed cannot."""
    name = record.get("species")
    return name == species.latin_name or species.latin_name in record.get("species_aliases", [])


def species_for_record(record):
    name = record.get("species")
    return Species._latin_names_to_species.get(name)


def records_for_selector(selector, species=None):
    records = catalogue()["dated"]
    matches = [r for r in records if _matches(r, selector)]
    if species is not None:
        matches = [r for r in matches if species_matches(species, r)]
    return matches


def infer_species(selector):
    matches = records_for_selector(selector)
    if not matches:
        # Historical numbered-only assembly names are valid assertions too.
        matches = [
            r for rows in catalogue()["numbered"].values() for r in rows.values()
            if _matches(r, selector)
        ]
    species_names = {r.get("species") for r in matches}
    if len(species_names) == 1:
        return species_for_record(matches[0])
    if len(species_names) > 1:
        raise ValueError("Assembly %r is ambiguous; specify species" % selector)
    return None


def validate_assembly(record, assembly=None, accession=None):
    if assembly is not None and not _matches(record, assembly):
        raise ValueError(
            "Assembly %r conflicts with selected %s (%s)"
            % (assembly, record.get("assembly_name") or record.get("assembly"),
               record.get("assembly_accession") or "accession unknown")
        )
    if accession is not None:
        check_accession(accession)
        actual = record.get("assembly_accession")
        if actual is None:
            raise ValueError("Cannot verify accession %s for this historical release" % accession)
        if actual != accession:
            raise ValueError("Assembly accession %s conflicts with selected %s" % (accession, actual))


def dated_assembly(date, species=None, assembly=None, accession=None, provider=None,
                   infer_historical=True):
    """Return a record and whether this is the registered species selection."""
    if isinstance(assembly, str):
        assembly = assembly.strip()
        if _ACCESSION.fullmatch(assembly.upper()):
            assembly = assembly.upper()
    if accession is not None:
        check_accession(accession)
    if provider is not None and (
        not isinstance(provider, str) or not re.fullmatch(r"[a-z][a-z0-9_]*", provider)
    ):
        raise ValueError("provider must be one lowercase directory name")
    selector = accession or assembly
    all_records = catalogue()["dated"]
    candidates = [r for r in all_records if _matches(r, selector)]
    if species is not None:
        compatible = [r for r in candidates if species_matches(species, r)]
        if selector and candidates and not compatible:
            raise ValueError("Assembly %r does not belong to %s" % (selector, species.latin_name))
        candidates = compatible
    if provider:
        candidates = [r for r in candidates if r["provider"] == provider]
    default = species.dated_releases if species is not None else None

    def default_match(r):
        return default and (r["assembly_accession"], r["provider"]) == tuple(default)

    dated = [r for r in candidates if date in r["dates"]]
    defaults = [r for r in candidates if default_match(r)]
    if assembly and not accession and not _ACCESSION.fullmatch(assembly) and not dated:
        available = [
            r for r in all_records if date in r["dates"]
            and species is not None and species_matches(species, r)
        ]
        if available:
            names = ", ".join(sorted({r["assembly"] for r in available}))
            raise ValueError("Dated release %s provides %s, not %s" % (date, names, assembly))
    if not selector and default and (
        not infer_historical or any(r in dated for r in defaults) or not dated
    ):
        candidates = defaults
    elif dated:
        candidates = dated
        preferred = [r for r in candidates if default_match(r)]
        if preferred:
            candidates = preferred
    # Several dates of one assembly/provider live in a single record.
    if len(candidates) > 1:
        choices = sorted({
            "%s (%s, %s)" % (r["assembly_name"], r["assembly_accession"], r["provider"])
            for r in candidates
        })
        raise ValueError("Ambiguous assembly/provider for %s; choose one of: %s"
                         % (date, "; ".join(choices)))
    if candidates:
        record = dict(candidates[0])
        validate_assembly(record, assembly, accession)
    elif accession or (isinstance(assembly, str) and _ACCESSION.fullmatch(assembly)):
        # Preserve arbitrary accession support. Unknown names/species cannot be
        # claimed as verified metadata; explicit API discovery can supply them.
        selected_accession = accession or assembly
        known = [r for r in all_records if r["assembly_accession"] == selected_accession]
        record = dict(known[0]) if known else dict(
            assembly=None, assembly_name=None, assembly_accession=selected_accession,
            species=species.latin_name if species else None, dates=[],
        )
        validate_assembly(record, assembly, selected_accession)
        record["provider"] = provider or "ensembl"
    elif selector:
        raise ValueError("Unknown assembly/provider selection %r; use a versioned accession or API discovery" % selector)
    elif default:
        known = [r for r in all_records if r["assembly_accession"] == default[0]]
        record = dict(known[0]) if known else dict(
            assembly=species.dated_reference_name, assembly_name=None,
            assembly_accession=default[0], provider=provider or default[1],
            species=species.latin_name, dates=[],
        )
        record["provider"] = provider or default[1]
    else:
        name = species.latin_name if species else "this selection"
        raise ValueError("No dated Ensembl releases for %s; specify an assembly accession" % name)
    if species is not None and record.get("species") and not species_matches(species, record):
        raise ValueError("Assembly does not belong to %s" % species.latin_name)
    standard = bool(default_match(record)) and accession is None and not (
        isinstance(assembly, str) and _ACCESSION.fullmatch(assembly)
    )
    return record, standard


def normalize_species(species, assembly=None, accession=None):
    if species is not None:
        return check_species_object(species)
    selector = accession or assembly
    return infer_species(selector) if selector else check_species_object("human")
