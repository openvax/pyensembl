"""Offline reference genetic-code assignments and their provenance."""

from copy import deepcopy
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
import json
import re


COMPARTMENTS = ("nuclear", "mitochondrial", "plastid")
_MITOCHONDRIAL = frozenset({"mt", "m", "chrmt", "chrm", "mito", "mtdna",
                          "mitochondrion", "mitochondrial", "mitochondrion_genome",
                          "mitochondrial_genome"})
_PLASTID = frozenset({"pt", "cp", "chrc", "chrpt", "chloroplast", "plastid", "apicoplast"})
_NUCLEAR = re.compile(r"(?:chr)?(?:[1-9][0-9]*|[1-4][lr]|[xywz]|i|ii|iii|iv|v|vi|vii|viii|ix|x|xi|xii|xiii|xiv|xv|xvi)\Z")


@lru_cache(maxsize=1)
def catalogue():
    return json.loads(files("pyensembl").joinpath("genetic_code_catalog.json").read_text())


def taxon_id(value):
    """Validate an NCBI identifier; unknown positive IDs remain unresolved."""
    if isinstance(value, str) and re.fullmatch(r"[1-9][0-9]*", value):
        value = int(value)
    if type(value) is not int or value < 1:
        raise ValueError("taxon_id must be a positive NCBI integer identifier")
    return value


def translation_table(value):
    """Validate an explicit NCBI table ID without interpreting its codons."""
    if type(value) is not int or str(value) not in catalogue()["code_names"]:
        raise ValueError("Unknown NCBI translation table: %r" % (value,))
    return value


def related_taxa(first, second):
    """Only equality or captured ancestor/descendant evidence is sufficient."""
    if first == second:
        return True
    taxa = catalogue()["taxa"]
    return (first in taxa.get(str(second), {}).get("lineage", [])
            or second in taxa.get(str(first), {}).get("lineage", []))


def _check_compartments(values):
    if not isinstance(values, dict) or any(
        not isinstance(contig, str) or not contig or compartment not in COMPARTMENTS
        for contig, compartment in values.items()
    ):
        raise ValueError("contig_compartments must map contig names to nuclear, mitochondrial or plastid")
    return dict(values)


def resolve_reference_metadata(value=None, *, species=None, reference_name=None,
                               source_taxon_id=None, species_taxon_id=None, taxon_source=None):
    """Pin a reference record; saved records do not refresh their assignments."""
    from .species import Species, find_species_by_name

    if value is not None and not isinstance(value, dict):
        raise TypeError("reference_metadata must be a dictionary")
    value = deepcopy(value or {})
    if value.get("schema_version") == 1:
        required = {"taxon_id", "species_taxon_id", "translation_tables",
                    "translation_table_sources", "contig_compartments", "taxonomy_source"}
        if not required.issubset(value):
            raise ValueError("Incomplete saved reference metadata")
        for key in ("taxon_id", "species_taxon_id"):
            if value.get(key) is not None:
                value[key] = taxon_id(value[key])
        value["contig_compartments"] = _check_compartments(value.get("contig_compartments", {}))
        tables = value.get("translation_tables", {})
        if not isinstance(tables, dict) or set(tables) != set(COMPARTMENTS):
            raise ValueError("Saved reference metadata must contain all three compartments")
        for code in tables.values():
            if code is not None:
                translation_table(code)
        sources = value["translation_table_sources"]
        if (not isinstance(sources, dict) or set(sources) != set(COMPARTMENTS)
                or any(source not in ("taxonomy", "reference", "unresolved") for source in sources.values())):
            raise ValueError("Invalid saved translation-table provenance")
        if not isinstance(value["taxonomy_source"], dict) or not isinstance(value["taxonomy_source"].get("url"), str):
            raise ValueError("Invalid saved taxonomy provenance")
        return value
    if "schema_version" in value:
        raise ValueError("Unsupported reference metadata schema")
    unknown = set(value) - {"taxon_id", "species_taxon_id", "translation_tables", "contig_compartments", "source_url"}
    if unknown:
        raise ValueError("Unknown reference metadata fields: %s" % ", ".join(sorted(unknown)))
    if species is None and reference_name is not None:
        species = Species._reference_names_to_species.get(reference_name)
    if isinstance(species, str):
        try:
            species = find_species_by_name(species)
        except ValueError:
            species = None
    data = catalogue()
    registered = data["species"].get(getattr(species, "latin_name", None), {})
    registered_id = registered.get("taxon_id")
    identifiers = [taxon_id(item) if item is not None else None for item in
                   (value.get("taxon_id"), source_taxon_id, registered_id,
                    value.get("species_taxon_id"), species_taxon_id)]
    selected = next((item for item in identifiers[:3] if item is not None), None)
    for other in identifiers:
        if other is not None and selected is not None and not related_taxa(selected, other):
            raise ValueError("Reference species and taxon_id disagree: %s versus %s" % (selected, other))
    taxon = data["taxa"].get(str(selected), {})
    parent_species = taxon.get("species_taxon_id") or identifiers[3] or identifiers[4]
    for declared in identifiers[3:]:
        if declared is not None and parent_species is not None and declared != parent_species:
            raise ValueError("species_taxon_id disagrees with the verified species ancestor: %s versus %s"
                             % (declared, parent_species))
    tables = dict.fromkeys(COMPARTMENTS)
    tables.update(taxon.get("translation_tables", {}))
    sources = {compartment: "taxonomy" if code is not None else "unresolved"
               for compartment, code in tables.items()}
    overrides = value.get("translation_tables", {})
    if not isinstance(overrides, dict) or set(overrides) - set(COMPARTMENTS):
        raise ValueError("translation_tables must map known compartments to NCBI table IDs")
    for compartment, code in overrides.items():
        tables[compartment] = translation_table(code) if code is not None else None
        sources[compartment] = "reference"
    source_url = value.get("source_url")
    if source_url is not None and not isinstance(source_url, str):
        raise ValueError("source_url must be a string")
    return dict(
        schema_version=1, taxon_id=selected, species_taxon_id=parent_species,
        scientific_name=taxon.get("scientific_name"),
        translation_tables=tables, translation_table_sources=sources,
        contig_compartments=_check_compartments(value.get("contig_compartments", {})),
        source_url=source_url,
        taxon_source=(source_url if identifiers[0] is not None else
                      taxon_source if source_taxon_id is not None else registered.get("source_url")),
        taxonomy_source=dict(data["sources"]["taxonomy.xml"], retrieved_at=data["retrieved_at"]),
    )


@dataclass(frozen=True)
class GeneticCodeAssignment:
    """An assignment or an explicit unresolved result; no sequence translation.

    The selected compartment defaults to nuclear when no compartment is known.
    It describes code selection, not verified physical localization.
    """

    table_id: int | None
    compartment: str | None
    source: str
    taxon_id: int | None
    source_url: str | None = None
    reason: str | None = None


def compartment_for_contig(contig, declared):
    if contig is None:
        return None
    contig = str(contig)
    if contig in declared:
        return declared[contig]
    matches = {compartment for name, compartment in declared.items() if name.lower() == contig.lower()}
    if len(matches) > 1:
        raise ValueError("Conflicting reference compartments for contig " + contig)
    if matches:
        return matches.pop()
    name = contig.lower()
    if name in _MITOCHONDRIAL:
        return "mitochondrial"
    if name in _PLASTID:
        return "plastid"
    return "nuclear" if _NUCLEAR.fullmatch(name) else None


def genetic_code(metadata, contig=None, *, compartment=None, source_table=None,
                 annotation_source=None):
    inferred = compartment_for_contig(contig, metadata["contig_compartments"])
    if compartment is not None:
        if compartment not in COMPARTMENTS:
            raise ValueError("Unknown cellular compartment: %r" % compartment)
        if inferred is not None and inferred != compartment:
            raise ValueError("Requested compartment conflicts with reference contig metadata")
    # Apply the species' main table only after validating known compartments;
    # an unfamiliar scaffold name must not conflict with an explicit request.
    compartment = compartment or inferred or "nuclear"
    identifier = metadata["taxon_id"]
    if source_table is not None:
        return GeneticCodeAssignment(translation_table(source_table), compartment,
                                     "annotation", identifier, annotation_source)
    code = metadata["translation_tables"][compartment]
    source = metadata["translation_table_sources"][compartment]
    url = metadata.get("source_url") if source == "reference" else metadata["taxonomy_source"]["url"]
    reason = None if code is not None else "No validated %s genetic-code assignment" % compartment
    return GeneticCodeAssignment(code, compartment, source, identifier, url, reason)
