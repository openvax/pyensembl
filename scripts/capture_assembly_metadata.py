#!/usr/bin/env python3
"""Capture bounded upstream assembly metadata without touching genome caches."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import re
from urllib.request import Request, urlopen
import zlib

from pyensembl.ensembl_url_templates import make_gtf_url
from pyensembl.species import Species


def get(url, limit=2 * 1024 * 1024, range_request=False):
    headers = {"Range": "bytes=0-%d" % (limit - 1)} if range_request else {}
    with urlopen(Request(url, headers=headers), timeout=30) as response:
        payload = response.read(limit if range_request else limit + 1)
    if len(payload) > limit:
        raise ValueError("Metadata response exceeds size limit: %s" % url)
    return payload


@lru_cache(maxsize=None)
def databases(url):
    return re.findall(r'href="([^"/]+)/"', get(url).decode())


def capture(root, pair):
    species, release = pair
    path = root / "headers" / ("%s-%s.json" % pair)
    if path.exists():
        return
    url = make_gtf_url(release, species)
    payload = get(url, 4096, range_request=True)
    text = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(payload).decode()
    headers = []
    for line in text.splitlines():
        if not line.startswith("#"):
            break
        headers.append(line)
    row = dict(species=species, release=release, url=url, headers=headers,
               prefix_sha256=hashlib.sha256(payload).hexdigest())
    if not any("genome-build-accession" in line for line in headers):
        directory = url.split("/gtf/")[0] + "/mysql/"
        try:
            matches = [name for name in databases(directory)
                       if name.startswith(species + "_core_")]
            if len(matches) != 1:
                raise ValueError("Ambiguous/missing core database: %r" % matches)
            meta_url = directory + matches[0] + "/meta.txt.gz"
            contents = gzip.decompress(get(meta_url)).decode()
            meta = {}
            for line in contents.splitlines():
                columns = line.split("\t")
                for i, key in enumerate(columns[:-1]):
                    if key.startswith("assembly.") or key == "species.taxonomy_id":
                        meta[key] = columns[i + 1]
            row.update(meta=meta, meta_url=meta_url)
        except (OSError, ValueError) as error:
            row["meta_error"] = str(error)
            print("%s %s: %s" % (species, release, error), flush=True)
    path.write_text(json.dumps(row, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--workers", type=int, default=4, choices=range(1, 9))
    args = parser.parse_args()
    (args.directory / "headers").mkdir(parents=True, exist_ok=True)
    platform = args.directory / "species.json"
    if not platform.exists():
        platform.write_bytes(get(
            "https://ftp.ebi.ac.uk/pub/ensemblorganisms/species.json", 64 * 1024 * 1024,
        ))
    pairs = list(Species.all_species_release_pairs())
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(lambda pair: capture(args.directory, pair), pairs))
    print("Captured metadata for %d species/release combinations" % len(pairs))
