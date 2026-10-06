# API reference

Find a class below. To find a query method by what it returns, use the
[method overview](genome.md#find-a-method). For installation and runnable
examples, start on the [home page](../index.md). The command-line interface has
a [separate reference](cli.md).

| Interface | Purpose | Complete reference |
| --- | --- | --- |
| `EnsemblRelease` | Select a species and numbered release | [Select data](datasets.md) |
| `EnsemblAnnotation` | Select an assembly and dated annotation | [Select data](datasets.md) |
| `Genome` | Shared annotation queries and custom data | [Query a genome](genome.md) |
| `Gene`, `Transcript`, `Exon`, `Protein`, `Locus` | Annotated features and intervals | [Features](features.md) |
| `GeneNameAliases` | Explicit alias source | [Aliases](aliases.md) |
| `DownloadCache`, `Database`, `SequenceData` | Local data and indexing | [Data and caches](data.md) |
| Reference and species helpers | Assembly/species normalization and selection | [Helpers](helpers.md) |

Coordinates in the annotation API are one-based and inclusive. The optional
[pyfaidx](https://github.com/mdshw5/pyfaidx) reader exposed by `Genome.fasta`
uses zero-based half-open slicing. Consult [reference DNA](../guides/reference-dna.md)
before mixing those interfaces. Installed annotation queries read local data;
setup methods explicitly download and index files.
