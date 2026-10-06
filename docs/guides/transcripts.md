# Get protein and transcript sequences

These examples go further than the [home page](../index.md#read-protein-and-transcript-sequences)
with the TP53 transcript TP53-201. They use the human GRCh38 / Ensembl release
93 data [installed on the home page](../index.md#install) and run in one Python
session.

## Protein sequence <a id="proteins"></a>

```python
from pyensembl import EnsemblRelease

data = EnsemblRelease(93, species="human")
transcript = data.transcript_by_id("ENST00000269305")
protein = transcript.protein_sequence
print(transcript.versioned_id, transcript.protein.versioned_id, len(protein))
```

```text
ENST00000269305.8 ENSP00000269305.4 393
```

`protein_sequence` is the annotated translation from Ensembl's peptide file.
It normally has no stop symbol, but a few peptides, mostly from polymorphic
pseudogenes, contain `*` stop symbols. `versioned_id` adds this
release's version to a stable ID; `transcript.protein_id` is the stable protein
ID without it.

When you start from a protein ID, look up its sequence, transcript or gene
directly:

```python
print(data.protein_sequence(transcript.protein_id) == protein)
print(data.transcript_by_protein_id("ENSP00000269305").name)
print(data.gene_by_protein_id("ENSP00000269305").name)
```

```text
True
TP53-201
TP53
```

Each lookup also accepts a [versioned ID](features.md#versioned-ids) such as
`ENSP00000269305.4`. A version that differs from this release's, such as
`ENSP00000269305.3`, raises `ValueError` instead of returning this release's
protein.

Noncoding transcripts have no protein, so their `protein`, `protein_id` and
`protein_sequence` are `None`:

```python
noncoding = data.transcript_by_id("ENST00000505014")
print(noncoding.name, noncoding.biotype, noncoding.protein_sequence)
```

```text
TP53-210 retained_intron None
```

## Coding sequence <a id="coding-sequence-and-utrs"></a>

```python
cds = transcript.coding_sequence
print(len(cds), cds[:9], cds[-3:])
print(transcript.complete)
```

```text
1182 ATGGAGGAG TGA
True
```

The coding sequence runs from the first base of the start codon through the
stop codon, so its 1,182 bases are 393 codons plus the stop. `complete` is
`True` when the transcript has annotated three-base start and stop codons and a
coding length divisible by three.

## Transcript sequence and UTRs

`transcript.sequence` is the full spliced cDNA, oriented 5′ to 3′. For a
transcript with annotated start and stop codons, it is the 5′ UTR, the coding
sequence and the 3′ UTR joined together:

```python
utr5 = transcript.five_prime_utr_sequence
utr3 = transcript.three_prime_utr_sequence
print(len(transcript.sequence), len(utr5), len(cds), len(utr3))
print(transcript.sequence == utr5 + cds + utr3)
```

```text
2579 190 1182 1207
True
```

To include introns or flanking sequence, read [genomic DNA](reference-dna.md)
from the same assembly.

## Incomplete transcripts

A `protein_coding` biotype does not guarantee a complete coding sequence.
Some transcripts are fragments without an annotated start or stop codon:

```python
fragment = data.transcript_by_id("ENST00000576024")
print(fragment.name, fragment.biotype, fragment.complete)
print(fragment.contains_start_codon, fragment.coding_sequence)
print(fragment.protein_sequence[:10])
```

```text
TP53-216 protein_coding False
False None
XSPQPKKKPL
```

Without an annotated start codon, `coding_sequence` and
`five_prime_utr_sequence` are `None`. The peptide file can still hold a partial
translation, so check `complete` before analyzing codons or reading frames.

## Exons and genomic coordinates <a id="transcript-span-and-exons"></a><a id="inspect-a-transcript"></a>

### Exons

```python
print(transcript.contig, transcript.start, transcript.end, transcript.strand)
for exon in transcript.exons[:2]:
    print(exon.id, exon.start, exon.end)
print(len(transcript.exons))
```

```text
17 7668402 7687538 -
ENSE00001146308 7687377 7687538
ENSE00002667911 7676521 7676622
11
```

The transcript's span includes its introns. `transcript.exons` is in
transcription order, so on the minus strand the first exon has the highest
coordinates. Every `start` is still the lower coordinate, and coordinates are
one-based and inclusive.

### Codon positions and coding intervals <a id="codon-positions-and-genomic-coding-intervals"></a>

```python
print(transcript.start_codon_positions)
print(transcript.stop_codon_positions)
print(transcript.coding_sequence_position_ranges[:2])
```

```text
[7676592, 7676593, 7676594]
[7669609, 7669610, 7669611]
[(7669609, 7669690), (7670609, 7670715)]
```

These are genomic positions in ascending order. On the minus strand,
translation starts at the highest position of the start codon, 7676594.
`coding_sequence_position_ranges` lists the coding part of each exon, including
the stop codon, which Ensembl annotates separately from the CDS.

### Map a genomic position to the transcript

`spliced_offset` converts a genomic position inside an exon to a zero-based
index into `transcript.sequence`:

```python
offset = transcript.spliced_offset(7676594)
print(offset, transcript.sequence[offset:offset + 3])
data.close()
```

```text
190 ATG
```

The first base of the start codon is at offset 190, immediately after the 5′
UTR. Positions in introns or outside the transcript raise `ValueError`. The
[Transcript reference](../reference/features.md#pyensembl.Transcript) lists
every attribute.
