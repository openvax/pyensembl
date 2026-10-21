#!/usr/bin/env python3
"""Build the offline assembly catalogue from captured Ensembl metadata.

Inputs are bounded GTF-header captures (including archive meta-table fields
where headers are absent), and the new platform's species.json. Keeping the
capture stage separate makes catalogue updates reviewable and reproducible.

    python scripts/build_assembly_catalog.py CAPTURE_DIRECTORY

See docs/dev/assembly-catalogue.md for acquisition and verification.
"""

import argparse
import json
from pathlib import Path
import re

from pyensembl.species import Species


def build(root):
    numbered = {}
    accessions = {}
    gaps = []
    for path in sorted((root / "headers").glob("*.json")):
        captured = json.loads(path.read_text())
        species = Species._latin_names_to_species[captured["species"]]
        release = captured["release"]
        header = dict(
            line[2:].split(" ", 1) for line in captured.get("headers", [])
            if line.startswith("#!") and " " in line
        )
        meta = captured.get("meta", {})
        accession = header.get("genome-build-accession") or meta.get("assembly.accession")
        if accession:
            accession = accession.rsplit(":", 1)[-1]
        if accession and not re.fullmatch(r"GC[AF]_\d{9}\.[1-9]\d*", accession):
            raise ValueError("Unexpected accession %r in %s" % (accession, path))
        record = dict(
            species=species.latin_name,
            assembly=header.get("genome-version") or meta.get("assembly.default")
                     or species.which_reference(release),
            assembly_name=header.get("genome-build") or meta.get("assembly.name"),
            assembly_accession=accession,
            reference_name=species.which_reference(release),
            source=captured.get("meta_url") if not header else captured["url"],
        )
        numbered.setdefault(species.latin_name, {})[str(release)] = record
        if accession:
            accessions.setdefault(accession, []).append(record)
        if not record["assembly_name"] or not accession:
            gaps.append(dict(species=species.latin_name, release=release,
                             missing=[key for key in ("assembly_name", "assembly_accession")
                                      if not record[key]]))
    expected = set(Species.all_species_release_pairs())
    observed = {(s, int(r)) for s, rows in numbered.items() for r in rows}
    if expected != observed:
        raise ValueError("Missing/extra numbered captures: %r" % (expected ^ observed))
    platform = json.loads((root / "species.json").read_text())
    default_accessions = {
        s.dated_releases[0]: s for s in Species._latin_names_to_species.values()
        if s.dated_releases
    }
    dated = []
    for organism in platform["species"].values():
        for accession, assembly in organism["assemblies"].items():
            if accession not in accessions and accession not in default_accessions:
                continue
            history = accessions.get(accession, [])
            species = default_accessions.get(accession)
            if species is None:
                species = Species._latin_names_to_species[history[0]["species"]]
            current = next(
                (r for r in reversed(history) if r["species"] == species.latin_name), None
            )
            for provider, dates in assembly["genebuild_providers"].items():
                record = dict(
                    species=species.latin_name,
                    assembly=assembly["name"], assembly_name=assembly["name"],
                    assembly_accession=accession, provider=provider,
                    dates=sorted(dates),
                    reference_name=species.dated_reference_name
                                   if accession in default_accessions
                                   else current["reference_name"] if current else assembly["name"],
                    source="https://ftp.ebi.ac.uk/pub/ensemblorganisms/species.json",
                )
                # GRC patch naming is verified in release metadata; other
                # species' version suffixes are meaningful and remain intact.
                builds = {r["assembly"] for r in history}
                if len(builds) == 1:
                    record["assembly"] = builds.pop()
                elif record["assembly_name"].startswith(("GRCh", "GRCm")):
                    record["assembly"] = record["assembly_name"].split(".p")[0]
                if species.latin_name.startswith("canis_lupus_familiaris"):
                    record["species_aliases"] = ["canis_familiaris"]
                elif species.latin_name.startswith("mus_musculus_"):
                    record["species_aliases"] = ["mus_musculus"]
                dated.append(record)
    # Retain the independently verified historical dataset even if a future
    # species.json omits it, without discarding any newly published dates.
    historical = dict(
        species="homo_sapiens", assembly="GRCh37", assembly_name="GRCh37.p13",
        assembly_accession="GCA_000001405.14", provider="ensembl", dates=["2013_09"],
        reference_name="GRCh37",
        source="https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/14/ensembl/2013_09/geneset/genes.gtf.gz",
    )
    historical["dates"] = sorted(set(historical["dates"]).union(
        date for record in dated for date in record["dates"]
        if (record["assembly_accession"], record["provider"]) ==
        (historical["assembly_accession"], historical["provider"])
    ))
    dated = [r for r in dated if (r["assembly_accession"], r["provider"]) !=
             (historical["assembly_accession"], historical["provider"])]
    dated.append(historical)
    # Keep every registered default, even if the bulk catalogue omits it.
    found = {(r["assembly_accession"], r["provider"]) for r in dated}
    for species in Species._latin_names_to_species.values():
        if not species.dated_releases or species.dated_releases in found:
            continue
        accession, provider = species.dated_releases
        candidates = accessions.get(accession, [])
        record = dict(candidates[-1]) if candidates else dict(
            species=species.latin_name, assembly=species.dated_reference_name,
            assembly_name=None, assembly_accession=accession,
            reference_name=species.dated_reference_name, source=None,
        )
        record.update(provider=provider, dates=[])
        dated.append(record)
    return dict(schema_version=1, platform_last_updated=platform["last_updated"],
                numbered=numbered, dated=dated, metadata_gaps=gaps)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("captures", type=Path)
    parser.add_argument("--output", type=Path,
                        default=Path("pyensembl/assembly_catalog.json"))
    args = parser.parse_args()
    result = build(args.captures)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print("%d numbered records, %d dated records, %d metadata gaps" % (
        sum(map(len, result["numbered"].values())), len(result["dated"]),
        len(result["metadata_gaps"]),
    ))
