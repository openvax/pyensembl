# Load custom genomes

Use `Genome` with local paths or remote URLs for matched
[GTF](https://en.wikipedia.org/wiki/Gene_transfer_format) annotations and
[FASTA](https://en.wikipedia.org/wiki/FASTA_format) sequences. GTF conventions
vary between providers; custom-data support depends on the supplied fields.
Matching contig names alone do not prove that the files use the same assembly.

## Index local annotations

This template requires your own GTF and FASTA files. Replace the paths before
running it; transcript and peptide record IDs must match the GTF identifiers.
Omit sequence sources when you only need annotation lookups.

```python
from pyensembl import Genome

with Genome(
    reference_name="GRCh38",
    annotation_name="my_genome_features",
    annotation_version="1",
    gtf_path_or_url="/data/my_genome_features.gtf",
    transcript_fasta_paths_or_urls=["/data/transcripts.fa"],
    protein_fasta_paths_or_urls=["/data/proteins.fa"],
    cache_directory_path="/data/pyensembl-indexes/my_genome_features-1",
) as data:
    data.index()
    gene_names = data.gene_names_at_locus(contig=6, position=29945884)
```

Local files are used in place by default, with indexes in the cache. Set
`copy_local_files_to_cache=True` for an independent cached import;
`decompress_on_download=True` also applies to copied sources. Remote sources
require `data.download()` before `data.index()`.

## Match files and identifiers

Keep annotation source, assembly, version and file provenance with your
analysis. A GTF should identify genes and transcripts using `gene_id` and
`transcript_id`; sequence lookup uses matching FASTA record IDs. Ensembl
version-suffix fallback supports matching a stable identifier to one versioned
record, but ambiguous versions are not chosen arbitrarily. Gene-only GTF files
can support gene queries without transcript identifiers.

GFF3 requires conversion to GTF before indexing. Unsupported or missing GTF
attributes can limit which queries work. [Data inspection](cache.md#inspect-data-without-installing)
checks files and index readiness, not biological validity or assembly identity.

To read genomic DNA, add `genome_fasta_path_or_url="/data/reference.fa"` using
the same assembly. See [local reference DNA](reference-dna.md#local-fasta-files)
for compression, read-only source files and contig validation. For a standard
new-platform geneset, use [Ensembl's dated annotations](../ensembl-platform.md).
