# Read genomic DNA

Reference DNA lets you read any genomic interval, including introns,
intergenic regions and flanking sequence. It is optional: a normal
installation does not download a whole genome. Human DNA is about 1 GB to
download and takes several GB of disk once decompressed.

## Quick start

```sh
# Annotation and DNA, or just the DNA:
pyensembl install --release 93 --species human --with-genome-fasta
pyensembl install --release 93 --species human --only-genome-fasta
```

```python
from pyensembl import EnsemblRelease

data = EnsemblRelease(93, species="human", genome_fasta=True)
data.download_genome_fasta()  # does nothing if the DNA is already installed
with data:
    bases = data.sequence("7", 117_480_000, 117_480_100)
    tp53 = data.gene_by_id("ENSG00000141510")  # TP53, on the minus strand
    tp53_dna = data.sequence(
        tp53.contig, tp53.start, tp53.end, strand=tp53.strand
    )
```

`genome_fasta=True` only chooses the DNA. Nothing is downloaded until you call
`download_genome_fasta()` or `download()`, or run `pyensembl install`.
Python objects use reference DNA only when constructed with `genome_fasta`:
a plain `EnsemblRelease(93)` does not pick up DNA installed by the CLI, and its
error message names the call that does.

## Reading sequences

`sequence(contig, start, end, mask="upper", *, strand="+")`:

- Coordinates are **one-based and inclusive**, like the rest of PyEnsembl,
  with `1 <= start <= end <= contig length`.
- Bases are read from the plus strand; `strand="-"` returns the reverse
  complement, so a gene or transcript reads 5' to 3'.
- Contigs can be named as in the FASTA or as PyEnsembl reports them
  (`gene.contig`). If the names differ only by a `chr` prefix, the error
  suggests the right one.
- Results are uppercase; `mask="raw"` keeps soft-masked repeats in lowercase.
- Absent contigs and invalid ranges raise `ValueError`. Missing DNA raises
  `MissingGenomeFastaError`, a `ValueError` whose message explains how to
  install or enable it. Reads never download anything.

Related attributes and methods:

- `fasta` is a [pyfaidx](https://github.com/mdshw5/pyfaidx) reader with
  zero-based, half-open slices (`fasta[contig][start - 1:end].seq`), as used by
  Varcode. It needs the FASTA's own contig names, and is `None` when DNA is not
  configured or not installed.
- `genome_fasta_path` is the uncompressed FASTA on disk, or `None`.
- `download()` and `index()` include DNA when it is configured.
  `index_genome_fasta()` builds the DNA index before the first query needs it.
- `close()` closes the reader. Readers already handed out stay usable after
  `clear_cache()`.
- Attached DNA does not affect equality: genes and transcripts from the same
  release compare equal with or without it.

## Choosing DNA

| | Python (`EnsemblRelease`) | CLI (`pyensembl install`) |
|---|---|---|
| Ensembl's DNA | `genome_fasta=True` | `--with-genome-fasta` or `--only-genome-fasta` |
| A local FASTA | `genome_fasta="/data/ref.fa.gz"` | `--genome-fasta-path /data/ref.fa.gz` |
| Coverage | `genome_fasta_type="primary_assembly"` | `--genome-fasta-type primary_assembly` |
| Masking | `genome_fasta_mask="soft"` | `--masked soft` |

The default is unmasked **toplevel** DNA, which covers the patch and haplotype
contigs in Ensembl annotations. `primary_assembly` has the chromosomes and
unplaced/unlocalized sequences but no patches or haplotypes, and some older
releases and species don't provide it. Masking is `none`, `soft` (repeats in
lowercase), or `hard` (repeats replaced with `N`). See Ensembl's
[DNA file definitions](https://ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/dna/README).

## Local FASTA files

```python
data = EnsemblRelease(
    93, species="human", genome_fasta="/data/my_reference.fa.gz"
)

# Custom annotations can attach DNA too, from a path or URL:
from pyensembl import Genome
custom = Genome("custom", "my_annotations",
                genome_fasta_path_or_url="/data/reference.fa")
```

Plain FASTA files are read in place; gzip and BGZF files are decompressed into
PyEnsembl's cache on first use. Indexes always live in the cache, so read-only
source directories work and your own files and indexes are never modified.
`index()` warns about annotation contigs that are missing from a local FASTA,
but matching contig names don't prove that the assembly matches.

## Disk space and shared files <a id="managing-disk-space"></a><a id="how-the-shared-dna-cache-works"></a>

Compatible releases share one copy of Ensembl DNA. To delete DNA, prune
unused copies or share a cache with other users, see
[free disk space](cache.md#free-disk-space) and
[how reference DNA is stored](cache.md#how-reference-dna-is-stored).

## Upgrading from 2.11.0

| 2.11.0 | 2.12.0 and later |
|---|---|
| `EnsemblRelease(81, download_genome_fasta=True)` | `EnsemblRelease(81, genome_fasta=True)` |
| `EnsemblRelease(81, genome_fasta_path="/data/ref.fa")` | `EnsemblRelease(81, genome_fasta="/data/ref.fa")` |

The old keywords still work but emit a `DeprecationWarning`. Objects pickled or
serialized by 2.11.0 still load. CLI flags are unchanged.
