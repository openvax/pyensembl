# Genome queries and setup

These methods apply to `EnsemblRelease`, `EnsemblAnnotation` and custom
[`Genome`](../guides/custom-genomes.md) objects. Find a method in the tables
below; the [complete reference](#pyensembl.Genome) follows them. Coordinates are
one-based and inclusive.

## Find a method

### Genes

| Method | Returns |
| --- | --- |
| [`gene_by_id`](#pyensembl.Genome.gene_by_id) | One `Gene` |
| [`genes_by_name`](#pyensembl.Genome.genes_by_name) | Every `Gene` with a name, optionally including [aliases](../guides/aliases.md) |
| [`genes_at_locus`](#pyensembl.Genome.genes_at_locus) | Genes overlapping a position or interval |
| [`nearest_gene`](#pyensembl.Genome.nearest_gene) | `(distance, Gene)` for the closest gene |
| [`genes`](#pyensembl.Genome.genes) | All genes, filtered by contig, strand or biotype |
| [`gene_by_protein_id`](#pyensembl.Genome.gene_by_protein_id) | The gene encoding a protein |
| [`merged_gene_intervals`](#pyensembl.Genome.merged_gene_intervals) | Non-overlapping gene intervals on a contig |

`gene_ids`, `gene_names`, `gene_ids_at_locus`, `gene_names_at_locus`,
`gene_ids_of_gene_name`, `gene_name_of_gene_id` and similar methods return IDs
or names instead of objects.

### Transcripts

| Method | Returns |
| --- | --- |
| [`transcript_by_id`](#pyensembl.Genome.transcript_by_id) | One `Transcript` |
| [`transcripts_by_name`](#pyensembl.Genome.transcripts_by_name) | Transcripts with a name such as TP53-201 |
| [`transcripts_at_locus`](#pyensembl.Genome.transcripts_at_locus) | Transcripts overlapping a position or interval |
| [`nearest_transcript`](#pyensembl.Genome.nearest_transcript) | `(distance, Transcript)` for the closest transcript |
| [`transcripts`](#pyensembl.Genome.transcripts) | All transcripts, filtered by contig, strand or biotype |
| [`transcript_by_protein_id`](#pyensembl.Genome.transcript_by_protein_id) | The transcript encoding a protein |

`transcript_ids`, `transcript_ids_of_gene_id`, `transcript_ids_at_locus` and
similar methods return IDs or names. A gene's transcripts are also available as
`gene.transcripts`.

### Exons

| Method | Returns |
| --- | --- |
| [`exon_by_id`](#pyensembl.Genome.exon_by_id) | One `Exon` |
| [`exons_at_locus`](#pyensembl.Genome.exons_at_locus) | Exons overlapping a position or interval |
| [`exons`](#pyensembl.Genome.exons) | All exons, filtered by contig or strand |

`exon_ids_of_transcript_id`, `exon_ids_of_gene_id` and similar methods return
IDs. Use `transcript.exons` for exons in transcription order.

### Proteins and sequences

| Method | Returns |
| --- | --- |
| [`protein_sequence`](#pyensembl.Genome.protein_sequence) | Amino acids for a protein ID |
| [`protein_ids`](#pyensembl.Genome.protein_ids) | Protein IDs, filtered by contig or strand |
| [`transcript_sequence`](#pyensembl.Genome.transcript_sequence) | Spliced cDNA for a transcript ID |

Transcript objects also provide `protein_sequence`, `coding_sequence`,
`sequence` and UTR sequences; see [protein and transcript sequences](../guides/transcripts.md).

### Reference DNA

These need [reference DNA](../guides/reference-dna.md), configured with
`genome_fasta`.

| Method | Returns |
| --- | --- |
| [`sequence`](#pyensembl.Genome.sequence) | Genomic DNA for an interval, optionally reverse-complemented |
| [`download_genome_fasta`](#pyensembl.Genome.download_genome_fasta) | Downloads only the DNA |
| [`index_genome_fasta`](#pyensembl.Genome.index_genome_fasta) | Builds the DNA index before the first query |
| [`fasta`](#pyensembl.Genome.fasta) | A zero-based pyfaidx reader, or `None` |
| [`genome_fasta_path`](#pyensembl.Genome.genome_fasta_path) | The uncompressed FASTA path, or `None` |

### Setup and data

| Method | Returns |
| --- | --- |
| [`download`](#pyensembl.Genome.download), [`index`](#pyensembl.Genome.index) | Downloads and indexes files; existing files are kept |
| [`installed`](#pyensembl.Genome.installed) | Whether queries are ready without setup |
| [`inspect_data`](#pyensembl.Genome.inspect_data) | Offline report of source files and indexes |
| [`contigs`](#pyensembl.Genome.contigs) | Contig names in the annotation |
| [`close`](#pyensembl.Genome.close) | Closes open index files; `with` does this automatically |

See [install and manage data](../guides/cache.md) for the command-line
equivalents.

::: pyensembl.Genome
