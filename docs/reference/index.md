# API reference

The reference is generated from the current Python interfaces and docstrings.
Use the [first lookup tutorial](../getting-started.md) for a runnable example.

| Interface | Purpose | Complete reference |
| --- | --- | --- |
| `EnsemblRelease` | Numbered species/release selection | [Dataset constructors](datasets.md) |
| `EnsemblAnnotation` | Explicit new-platform assembly/date selection | [Dataset constructors](datasets.md) |
| `Genome` | Shared annotation queries and custom data | [Genome](genome.md) |
| `Gene`, `Transcript`, `Exon`, `Protein`, `Locus` | Annotated features and intervals | [Features](features.md) |
| `GeneNameAliases` | Explicit alias source | [Aliases](aliases.md) |
| `DownloadCache`, `Database`, `SequenceData` | Local data and indexing | [Data and caches](data.md) |
| Reference and species helpers | Assembly/species normalization and selection | [Helpers](helpers.md) |

Coordinates in the annotation API are one-based and inclusive. The optional
[pyfaidx](https://github.com/mdshw5/pyfaidx) reader exposed by `Genome.fasta`
uses zero-based half-open slicing. Consult [reference DNA](../guides/reference-dna.md)
before mixing those interfaces. Installed annotation queries read local data;
setup methods explicitly download and index files.

The [lookup overview](lookups.md) preserves the original method descriptions
and source links. The [CLI reference](cli.md) lists commands and setup options.
