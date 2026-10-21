#!/usr/bin/env python3
"""Explicitly capture bounded Ensembl/NCBI metadata into a new directory.

No genome downloads, cache changes or network access during catalogue builds.
"""

import argparse
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

from build_genetic_code_catalog import SOURCES, species_mapping


def get(url):
    with urlopen(Request(url, headers={"User-Agent": "PyEnsembl-metadata-capture"}), timeout=60) as response:
        payload = response.read(16 * 1024 * 1024 + 1)
    if len(payload) > 16 * 1024 * 1024:
        raise ValueError("Metadata response exceeds 16 MiB: " + url)
    return payload


def taxonomy_url(identifiers):
    return "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode(
        dict(db="taxonomy", id=",".join(map(str, sorted(identifiers))), retmode="xml")
    )


def capture(root):
    root.mkdir(parents=True, exist_ok=False)
    for name, url in SOURCES.items():
        (root / (name + ".gz")).write_bytes(gzip.compress(get(url), mtime=0))
    identifiers = {row["taxon_id"] for row in species_mapping(root).values()}
    first = get(taxonomy_url(identifiers))
    for taxon in ET.fromstring(first).findall("Taxon"):
        for ancestor in taxon.findall("LineageEx/Taxon"):
            if ancestor.findtext("Rank") == "species":
                identifiers.add(int(ancestor.findtext("TaxId")))
    url = taxonomy_url(identifiers)
    payload = get(url)
    returned = {int(node.findtext("TaxId")) for node in ET.fromstring(payload).findall("Taxon")}
    if returned != identifiers:
        raise ValueError("NCBI did not return all requested taxa: %r" % (identifiers - returned))
    (root / "taxonomy.xml.gz").write_bytes(gzip.compress(payload, mtime=0))
    (root / "capture.json").write_text(json.dumps(dict(
        retrieved_at=datetime.now(timezone.utc).isoformat(), taxonomy_url=url,
    ), indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    capture(parser.parse_args().directory)
