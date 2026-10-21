"""Explicit discovery through the new Ensembl metadata API.

This module is never called by ordinary construction or deserialization.
"""

import json
import re
from urllib.parse import urlsplit
from urllib.request import Request, urlopen
from uuid import UUID

from .assembly import check_accession, records_for_selector
from .ensembl_versions import normalize_release
from .species import find_species_by_name
from .genetic_code import taxon_id


METADATA_API = "https://www.ensembl.org/api/metadata"
_DATASET_PATH = re.compile(
    r"/pub/ensemblorganisms/(GC[AF])/(\d{3})/(\d{3})/(\d{3})/([1-9]\d*)/"
    r"([a-z][a-z0-9_]*)/(\d{4}_\d{2})/(geneset|genome)/?\Z"
)


def _get_json(path):
    request = Request(METADATA_API + path, headers={"Accept": "application/json"})
    with urlopen(request, timeout=30) as response:
        # Metadata should be small. A proxy HTML page or unbounded response is
        # not a dataset; do not download annotation files through this client.
        payload = response.read(2 * 1024 * 1024 + 1)
    if len(payload) > 2 * 1024 * 1024:
        raise ValueError("Ensembl metadata response exceeds 2 MiB")
    try:
        return json.loads(payload)
    except (ValueError, UnicodeError) as error:
        raise ValueError("Invalid JSON from Ensembl metadata API") from error


def _dataset(url, kind):
    if not isinstance(url, str):
        raise ValueError("Missing Ensembl dataset URL")
    parsed = urlsplit(url)
    match = _DATASET_PATH.fullmatch(parsed.path)
    if (parsed.scheme != "https" or parsed.netloc != "ftp.ebi.ac.uk"
            or parsed.query or parsed.fragment or match is None
            or match.group(8) != kind):
        raise ValueError("Unexpected Ensembl %s URL: %s" % (kind, url))
    prefix, a, b, c, version, provider, date, _ = match.groups()
    return "%s_%s%s%s.%s" % (prefix, a, b, c, version), provider, normalize_release(date)


def resolve_genome(identifier):
    if not isinstance(identifier, str):
        raise TypeError("Expected an Ensembl genome UUID or versioned assembly accession")
    accession = None
    if identifier.startswith(("GCA_", "GCF_")):
        accession = check_accession(identifier)
        explained = _get_json("/genome/%s/explain" % accession)
        if not isinstance(explained, dict) or explained.get("is_suppressed"):
            raise ValueError("Ensembl assembly is unavailable or suppressed")
        if not isinstance(explained.get("assembly"), dict) or (
            explained["assembly"].get("accession_id") != accession
        ):
            raise ValueError("Ensembl API returned a different assembly accession")
        identifier = explained.get("genome_id")
    try:
        genome_id = str(UUID(identifier))
    except (ValueError, TypeError, AttributeError) as error:
        raise ValueError("Expected a valid Ensembl genome UUID") from error
    details = _get_json("/genome/%s/details" % genome_id)
    links = _get_json("/genome/%s/ftplinks" % genome_id)
    if not isinstance(details, dict) or details.get("genome_id") != genome_id:
        raise ValueError("Ensembl API returned a different genome UUID")
    if details.get("is_suppressed"):
        raise ValueError("Ensembl genome is suppressed")
    assembly = details.get("assembly", {})
    if not isinstance(assembly, dict):
        raise ValueError("Ensembl metadata has no assembly record")
    actual = check_accession(assembly.get("accession_id"))
    if accession is not None and actual != accession:
        raise ValueError("Ensembl API returned a different assembly accession")
    name = assembly.get("name")
    if not isinstance(name, str) or not name:
        raise ValueError("Ensembl metadata has no assembly name")
    if not isinstance(links, list) or any(not isinstance(row, dict) for row in links):
        raise ValueError("Ensembl metadata has no dataset links")
    genesets = [row.get("url") for row in links if row.get("dataset") == "genebuild"]
    if len(genesets) != 1:
        raise ValueError("Expected exactly one Ensembl geneset; specify a genome UUID")
    accession, provider, date = _dataset(genesets[0], "geneset")
    if accession != actual:
        raise ValueError("Ensembl geneset URL conflicts with assembly metadata")
    for row in links:
        if row.get("dataset") == "assembly":
            if _dataset(row.get("url"), "genome") != (accession, provider, date):
                raise ValueError("Ensembl genome and geneset URLs select different datasets")
    annotation_date = details.get("annotation_date")
    if annotation_date is not None and not isinstance(annotation_date, str):
        raise ValueError("Ensembl metadata has an invalid annotation date")
    if annotation_date and annotation_date.replace("-", "_") != date:
        raise ValueError("Ensembl annotation date conflicts with its geneset URL")
    known = records_for_selector(accession)
    record = dict(known[0]) if known else {}
    if not record.get("species") and details.get("scientific_name"):
        try:
            record["species"] = find_species_by_name(details["scientific_name"]).latin_name
        except ValueError:
            pass
    # The exact official name is supplied by this selected API record. For an
    # unregistered assembly it is also the best known build identifier.
    record.update(
        assembly=record.get("assembly") or name,
        assembly_name=name, assembly_accession=accession, provider=provider,
        genome_id=genome_id, website_release=details.get("release", {}).get("name"),
        scientific_name=details.get("scientific_name"),
        source=METADATA_API + "/genome/%s/details" % genome_id, dates=[date],
    )
    for incoming, field in (("taxonomy_id", "taxon_id"), ("species_taxonomy_id", "species_taxon_id")):
        if details.get(incoming) is not None:
            record[field] = taxon_id(details[incoming])
    return record, date
