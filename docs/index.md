# PyEnsembl

PyEnsembl queries [Ensembl](https://www.ensembl.org/) gene annotations and
sequences from local indexes. Select an assembly and annotation, install its
matched [GTF](https://en.wikipedia.org/wiki/Gene_transfer_format) and
[FASTA](https://en.wikipedia.org/wiki/FASTA_format) files, then look up genes,
transcripts, exons and proteins in Python.

Start with [installation and a first lookup](getting-started.md). The tutorial
pins human GRCh38 to Ensembl release 93 and explains the returned identifiers,
coordinates, strand and sequence lengths.

For other data, [choose an assembly](guides/assembly-selection.md), use the
[new Ensembl platform](ensembl-platform.md), or load
[custom genomes](guides/custom-genomes.md). The [API reference](reference/index.md)
contains complete interfaces; task guides explain aliases, cache readiness and
optional reference DNA.

Gene names can match several loci, and annotation versions can change
coordinates and sequences. Keep the chosen data version with your analysis.
Reference DNA is optional and must match the annotation's assembly.
