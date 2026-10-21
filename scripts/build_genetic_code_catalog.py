#!/usr/bin/env python3
"""Build a deterministic offline genetic-code catalogue from captured metadata.

Inputs: Ensembl species TSVs, taxonomy.xml, gc.prt and capture.json. Each input
may instead be gzip compressed. Capture is explicit and separate from builds.
"""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from pyensembl.assembly import numbered_assembly
from pyensembl.species import Species


SOURCES = {
    "vertebrates.tsv": "https://ftp.ensembl.org/pub/release-116/species_EnsemblVertebrates.txt",
    **{name + ".tsv": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-63/%s/species_Ensembl%s.txt"
       % (name, name.title()) for name in ("plants", "fungi", "protists", "metazoa")},
    "gc.prt": "https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt",
}


def read(root, name):
    path = root / name
    return path.read_bytes() if path.exists() else gzip.decompress(path.with_suffix(path.suffix + ".gz").read_bytes())


def species_mapping(root):
    rows = {name: list(csv.DictReader(read(root, name + ".tsv").decode().splitlines(), delimiter="\t"))
            for name in ("vertebrates", "plants", "fungi", "protists", "metazoa")}
    result = {}
    for species in Species._latin_names_to_species.values():
        release = max(end for _, end in species.reference_assemblies.values())
        division = species.division if species.ensembl_genomes else "vertebrates"
        matches = [row for row in rows[division] if row["species"] == species.ensembl_name(release)]
        if not matches:
            accession = numbered_assembly(species, release)["assembly_accession"]
            matches = [row for row in rows[division] if accession and row["assembly_accession"] == accession]
        if len(matches) != 1:
            raise ValueError("Expected one Ensembl taxon mapping for " + species.latin_name)
        row = matches[0]
        result[species.latin_name] = dict(taxon_id=int(row["taxonomy_id"]),
                                        source_url=SOURCES[division + ".tsv"])
    return result


def build(root):
    root = Path(root)
    capture = json.loads(read(root, "capture.json"))
    code_names = {}
    for block in re.findall(r"\{([^{}]+)\}", read(root, "gc.prt").decode()):
        identifier = re.search(r"\bid\s+(\d+)", block)
        if identifier:
            code_names[identifier[1]] = re.search(r'name\s+"([^"]+)"', block)[1]

    def code(value):
        if value in (None, "0"):
            return None
        if value not in code_names:
            raise ValueError("Unknown NCBI translation table: " + value)
        return int(value)

    taxa = {}
    for taxon in ET.fromstring(read(root, "taxonomy.xml")).findall("Taxon"):
        identifier = int(taxon.findtext("TaxId"))
        lineage = taxon.findall("LineageEx/Taxon")
        species = identifier if taxon.findtext("Rank") == "species" else next(
            (int(node.findtext("TaxId")) for node in reversed(lineage) if node.findtext("Rank") == "species"), None)
        plastid = next((node.findtext("PropValueInt") for node in taxon.findall("Properties/Property")
                        if node.findtext("PropName") == "pgcode"), None)
        aliases = {taxon.findtext("ScientificName")}
        for tag in ("Synonym", "EquivalentName", "GenbankCommonName", "CommonName"):
            aliases.update(node.text for node in taxon.findall("OtherNames/" + tag))
        taxa[str(identifier)] = dict(
            scientific_name=taxon.findtext("ScientificName"),
            aliases=sorted(name.lower().replace(" ", "_") for name in aliases if name),
            parent_id=int(taxon.findtext("ParentTaxId")), species_taxon_id=species,
            lineage=[int(node.findtext("TaxId")) for node in lineage],
            translation_tables=dict(nuclear=code(taxon.findtext("GeneticCode/GCId")),
                                    mitochondrial=code(taxon.findtext("MitoGeneticCode/MGCId")),
                                    plastid=code(plastid)),
        )
    species = species_mapping(root)
    for name, row in species.items():
        if str(row["taxon_id"]) not in taxa:
            raise ValueError("Missing captured NCBI taxon for " + name)
    source_urls = dict(SOURCES, **{"taxonomy.xml": capture["taxonomy_url"]})
    sources = {name: dict(url=url, sha256=hashlib.sha256(read(root, name)).hexdigest())
               for name, url in source_urls.items()}
    return dict(schema_version=1, retrieved_at=capture["retrieved_at"],
                sources=sources, code_names=code_names, species=species, taxa=taxa)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("captures", type=Path)
    parser.add_argument("--output", type=Path, default=Path("pyensembl/genetic_code_catalog.json"))
    args = parser.parse_args()
    result = build(args.captures)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print("%d registered species, %d NCBI taxa" % (len(result["species"]), len(result["taxa"])))
