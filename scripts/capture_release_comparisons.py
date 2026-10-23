"""Reproduce scoped numbered/dated GTF comparisons from complete source files.

Run: python scripts/capture_release_comparisons.py --output /tmp/ensembl-evidence
Downloads up to 1 GiB per source into that explicit directory (about 1.2 GiB
for the current 20 sources). Never reads or writes an installed genome cache.
Results compare GTF annotations, not cDNA, peptide or DNA sequence contents.
"""

import argparse
import concurrent.futures
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
from urllib.request import urlopen

from pyensembl import EnsemblRelease

ROOT = None
ATTR = re.compile(r'(\w+)\s+"([^"]*)"')
PAIRS = [
    ('human', 110, '2023_03'), ('human', 114, '2024_11'),
    ('human', 115, '2025_12'), ('human', 116, '2026_04'),
    ('human', 75, '2013_09'), ('mouse', 116, '2026_04'),
    ('dog', 116, '2020_10'), ('yeast', 116, '2018_10'),
    ('fly', 116, '2022_07'), ('arabidopsis_thaliana', 63, '2010_09'),
]


def acquire(specification):
    species, release = specification
    genome = EnsemblRelease(release, species=species)
    url = genome.gtf_url
    path = ROOT / ('%s-%s.gtf.gz' % (species, release))
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    if not path.exists():
        temporary = path.with_suffix('.part')
        with urlopen(url, timeout=90) as response, temporary.open('wb') as output:
            count = 0
            while block := response.read(1024 * 1024):
                count += len(block)
                if count > 1024**3:
                    raise ValueError('GTF exceeds 1 GiB download bound: ' + url)
                output.write(block)
        temporary.replace(path)
    with path.open('rb') as source:
        while block := source.read(4 * 1024 * 1024):
            md5.update(block)
            sha256.update(block)
    result = dict(species=species, release=release, source=url,
                  assembly_accession=genome.assembly_accession,
                  path=str(path), bytes=path.stat().st_size,
                  md5=md5.hexdigest(), sha256=sha256.hexdigest())
    if isinstance(release, str):
        manifest_url = url.rsplit('/', 1)[0] + '/md5sum.txt'
        with urlopen(manifest_url, timeout=30) as response:
            manifest = response.read(2*1024*1024)
        matches = [line.split()[0] for line in manifest.decode().splitlines()
                   if len(line.split()) == 2 and line.split()[1].lstrip('*') == url.rsplit('/', 1)[1]]
        assert matches == [result['md5']], (url, matches, result['md5'])
        (ROOT / ('%s-%s-md5sum.txt' % (species, release))).write_bytes(manifest)
        result.update(manifest_url=manifest_url, manifest_sha256=hashlib.sha256(manifest).hexdigest())
    (ROOT / ('%s-%s-acquired.json' % (species, release))).write_text(json.dumps(result, indent=2))
    print('Acquired', species, release, result['bytes'], flush=True)
    return result


def summarize(record):
    global ROOT
    ROOT = Path(record["path"]).parent
    genes, transcripts, rows = set(), set(), []
    headers, counts, sample = [], {}, []
    with gzip.open(record['path'], 'rt') as source:
        for line in source:
            if line.startswith('#'):
                headers.append(line.rstrip())
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError('Not GTF: ' + record['path'])
            pairs = ATTR.findall(fields[8])
            if ATTR.sub("", fields[8]).replace(";", "").strip():
                raise ValueError("Unsupported GTF attribute syntax: " + fields[8])
            attributes = dict(pairs)
            kind = fields[2]
            counts[kind] = counts.get(kind, 0) + 1
            if kind == 'gene':
                genes.add((attributes['gene_id'], attributes.get('gene_version', '')))
            if kind == 'transcript':
                transcripts.add((attributes['transcript_id'], attributes.get('transcript_version', '')))
            normalized = json.dumps([*fields[:8], sorted(pairs)], separators=(',', ':'))
            rows.append(hashlib.sha256(normalized.encode()).digest())
            if len(sample) < 100:
                sample.append(line)
    def digest(values):
        return hashlib.sha256(json.dumps(sorted(values), separators=(',', ':')).encode()).hexdigest()
    rows.sort()
    result = {k:v for k,v in record.items() if k != 'path'}
    result.update(genes=len(genes), transcripts=len(transcripts),
                  gene_ids_sha256=digest(genes), transcript_ids_sha256=digest(transcripts),
                  normalized_gtf_sha256=hashlib.sha256(b''.join(rows)).hexdigest(),
                  feature_counts=counts, headers=headers,
                  retrieved_at=datetime.now(timezone.utc).isoformat())
    stem='%s-%s' % (record['species'], record['release'])
    for kind, identifiers in [('genes', genes), ('transcripts', transcripts)]:
        (ROOT / (stem+'-'+kind+'.json.gz')).write_bytes(
            gzip.compress(json.dumps(sorted(identifiers), separators=(',', ':')).encode(), mtime=0))
    (ROOT / (stem+'-summary.json')).write_text(json.dumps(result, indent=2))
    (ROOT / (stem+'-sample.gtf')).write_text('\n'.join(headers)+'\n'+''.join(sample))
    print('Summarized', stem, len(genes), len(transcripts), flush=True)
    return result


def compare(summaries):
    by_key = {(r["species"], str(r["release"])): r for r in summaries}
    comparisons = []
    for species, numbered, dated in PAIRS:
        left, right = by_key[species, str(numbered)], by_key[species, dated]
        def identifiers(release, kind):
            path = ROOT / ("%s-%s-%s.json.gz" % (species, release, kind))
            return set(map(tuple, json.loads(gzip.decompress(path.read_bytes()))))
        genes = [identifiers(r, "genes") for r in (numbered, dated)]
        transcripts = [identifiers(r, "transcripts") for r in (numbered, dated)]
        same_ids = genes[0] == genes[1] and transcripts[0] == transcripts[1]
        same_rows = left["normalized_gtf_sha256"] == right["normalized_gtf_sha256"]
        subset = genes[1] <= genes[0] and transcripts[1] < transcripts[0]
        relationship = ("same-gtf-records" if same_rows else "same-identifiers" if same_ids
                        else "transcript-subset" if subset else "different")
        comparison = dict(
            species=left["species"], assembly_accession=right["assembly_accession"],
            provider=EnsemblRelease(dated, species=species).provider,
            include_alt="genes-including_alt" in right["source"],
            numbered=numbered, dated=dated, relationship=relationship,
            numbered_source=left, dated_source=right,
        )
        for kind, sets in (("genes", genes), ("transcripts", transcripts)):
            bare = [{identifier for identifier, version in values} for values in sets]
            comparison[kind] = dict(
                only_numbered=len(sets[0] - sets[1]), only_dated=len(sets[1] - sets[0]),
                shared=len(sets[0] & sets[1]),
                unversioned_only_numbered=len(bare[0] - bare[1]),
                unversioned_only_dated=len(bare[1] - bare[0]),
                unversioned_shared=len(bare[0] & bare[1]),
                numbered_missing_version=sum(not version for _, version in sets[0]),
                dated_missing_version=sum(not version for _, version in sets[1]),
            )
        comparisons.append(comparison)
    return dict(
        schema_version=1,
        method="Quoted GTF attributes; gene/transcript ID and version pairs; SHA256 of sorted row hashes including columns 1-8 and sorted attribute pairs. Headers and row order excluded; duplicates retained. No FASTA equivalence claim.",
        comparisons=comparisons, preferred=[],
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summaries-only", action="store_true", help="Compare existing complete captures without downloading or reparsing")
    args = parser.parse_args()
    ROOT = args.output
    ROOT.mkdir(parents=True, exist_ok=True)
    selected = list(dict.fromkeys((s, r) for s, left, right in PAIRS for r in (left, right)))
    if args.summaries_only:
        summaries = [json.loads((ROOT / ("%s-%s-summary.json" % pair)).read_text()) for pair in selected]
    else:
        with concurrent.futures.ThreadPoolExecutor(4) as pool:
            acquired = list(pool.map(acquire, selected))
        with concurrent.futures.ProcessPoolExecutor(2) as pool:
            summaries = list(pool.map(summarize, acquired))
    (ROOT / "release_comparisons.json").write_text(json.dumps(compare(summaries), indent=2) + "\n")
