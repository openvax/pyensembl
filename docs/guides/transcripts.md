# Exons and coding sequences

A transcript's exons, coding sequence and UTRs are available from its
annotation. These examples use the human GRCh38 / Ensembl release 93 data
installed on the [home page](../index.md#install) and the TP53 transcript
TP53-201. They run in one Python session.

## Transcript span and exons <a id="inspect-a-transcript"></a>

```python
from pyensembl import EnsemblRelease

data = EnsemblRelease(93, species="human")
transcript = data.transcript_by_id("ENST00000269305")
print(transcript.id, transcript.version, transcript.name, transcript.biotype)
print(transcript.contig, transcript.start, transcript.end, transcript.strand)
for exon in transcript.exons[:2]:
    print(exon.id, exon.start, exon.end)
print(len(transcript.exons))
```

```text
ENST00000269305 8 TP53-201 protein_coding
17 7668402 7687538 -
ENSE00001146308 7687377 7687538
ENSE00002667911 7676521 7676622
11
```

The transcript's span includes its introns. `transcript.exons` is in
transcription order, so on the minus strand the first exon has the highest
coordinates. Every `start` is still the lower coordinate, and coordinates are
one-based and inclusive.

## Coding sequence and UTRs

`transcript.sequence` is the spliced cDNA. It divides into the 5′ UTR, the
coding sequence and the 3′ UTR:

```python
utr5 = transcript.five_prime_utr_sequence
cds = transcript.coding_sequence
utr3 = transcript.three_prime_utr_sequence
print(len(utr5), len(cds), len(utr3), len(transcript.sequence))
print(cds[:9], cds[-3:], len(transcript.protein_sequence))
print(transcript.complete)
```

```text
190 1182 1207 2579
ATGGAGGAG TGA 393
True
```

The coding sequence runs from the first base of the start codon through the
stop codon, so its 1,182 bases encode 393 amino acids plus the stop.
`complete` is `True` when the transcript has annotated start and stop codons
and a coding length divisible by three.

## Codon positions and genomic coding intervals

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

## Map a genomic position to the transcript

`spliced_offset` converts a genomic position inside an exon to a zero-based
index into `transcript.sequence`:

```python
offset = transcript.spliced_offset(7676594)
print(offset, transcript.sequence[offset:offset + 3])
```

```text
190 ATG
```

The first base of the start codon is at offset 190, immediately after the 5′
UTR. Positions in introns or outside the transcript raise `ValueError`.

## Proteins

```python
print(transcript.protein_id)
print(data.transcript_by_protein_id("ENSP00000269305").name)
print(data.gene_by_protein_id("ENSP00000269305").name)
```

```text
ENSP00000269305
TP53-201
TP53
```

`transcript.protein_sequence` is the annotated translation from Ensembl's
peptide file. Noncoding transcripts have no `protein_id` and return `None`.

## Incomplete transcripts

A `protein_coding` biotype does not guarantee a complete coding sequence.
Some transcripts are fragments without an annotated start or stop codon:

```python
fragment = data.transcript_by_id("ENST00000576024")
print(fragment.name, fragment.biotype, fragment.complete)
print(fragment.contains_start_codon, fragment.coding_sequence)
print(fragment.protein_sequence[:10])
data.close()
```

```text
TP53-216 protein_coding False
False None
XSPQPKKKPL
```

Without an annotated start codon, `coding_sequence` and
`five_prime_utr_sequence` are `None`. The peptide file can still hold a partial
translation, so check `complete` before analyzing codons or reading frames.
The [Transcript reference](../reference/features.md#pyensembl.Transcript)
lists every attribute. To include introns or flanking sequence, read
[genomic DNA](reference-dna.md) from the same assembly.
