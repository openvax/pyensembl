# Command-line reference

Run `pyensembl ACTION [options]`. [Installation](../index.md#install),
[assembly selection](../guides/assembly-selection.md), [data inspection](../guides/cache.md)
and [reference DNA](../guides/reference-dna.md) explain the main workflows.
Dated releases of the new platform use `--release YYYY_MM`, e.g.
`--release 2026_04`; other new-platform datasets use the Python constructor
described in the [platform guide](../ensembl-platform.md).

## Commands and options

This is the complete help output for the documented package version. Run
`pyensembl --help` for the version installed in your environment.

```text
usage: pyensembl ACTION [options]

Install and manage the genome data PyEnsembl uses.

positional arguments:
  ACTION                one of the actions listed below

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -v, --verbose         Show detailed progress, including download and
                        database steps
  --json                With inspect, print a JSON inventory instead of a file
                        table
  --overwrite           Force download and indexing even if files already
                        exist locally
  --reference-name REFERENCE_NAME
                        Reference assembly, e.g. GRCh37 (case-insensitive for
                        Ensembl). Selects its newest supported release unless
                        --release is given; with custom source files, names
                        the custom reference.

Ensembl release options:
  --release RELEASE [RELEASE ...]
                        Ensembl release version(s), numbered or a YYYY_MM
                        annotation date on the new Ensembl platform; required
                        for deletion (install default=newest supported release
                        for --reference-name, otherwise for each species)
  --species SPECIES [SPECIES ...]
                        Which species to download Ensembl data for
                        (default=inferred from --reference-name, otherwise
                        human)
  --custom-mirror CUSTOM_MIRROR
                        URL and directory to use instead of the default
                        Ensembl FTP server

Custom genome options:
  --annotation-name ANNOTATION_NAME
                        Name of annotation source (e.g. refseq)
  --annotation-version ANNOTATION_VERSION
                        Version of annotation database
  --gtf GTF             URL or local path to a GTF file containing
                        annotations.
  --transcript-fasta TRANSCRIPT_FASTA
                        URL or local path to a FASTA files containing the
                        transcript data. This option can be specified multiple
                        times for multiple FASTA files.
  --protein-fasta PROTEIN_FASTA
                        URL or local path to a FASTA file containing protein
                        data.
  --shared-prefix SHARED_PREFIX
                        Add this prefix to URLs or paths specified by --gtf,
                        --transcript-fasta, --protein-fasta

Optional reference DNA:
  --with-genome-fasta   Also download and index reference DNA (several GB for
                        human)
  --only-genome-fasta   Download/index only reference DNA, without annotation
                        or transcript data
  --genome-fasta-path GENOME_FASTA_PATH
                        Attach a local reference FASTA; custom Genome sources
                        also accept a URL
  --genome-fasta-type {toplevel,primary_assembly}
                        Toplevel includes patch/haplotype contigs (default)
  --masked {none,soft,hard}
                        Masking of downloaded reference DNA (default: none)
  --check-genome-fasta  With list or inspect, check DNA indexes without
                        downloading or rebuilding
  --dry-run             With prune, report candidates without deleting files

actions:
  install             download and index genome data (skips what is already done)
  list                show installed genomes, whether they are indexed, and their DNA
  inspect             inspect selected source files, indexes, and download provenance offline
  available           show supported species, assemblies, and Ensembl releases
  delete-index-files  delete indexes, keeping downloaded files (needs --release)
  delete-all-files    delete all of a genome's local data (needs --release)
  prune               delete shared reference DNA that no installed release uses

examples:
  pyensembl install --release 75 77                     human releases 75 and 77
  pyensembl install --release 110 --species mouse       a mouse release
  pyensembl install --release 2026_04                   a dated release (new Ensembl platform)
  pyensembl install --reference-name GRCh37             newest release for GRCh37
  pyensembl install --release 110 --with-genome-fasta   also install reference DNA
  pyensembl install --reference-name GRCh38 --annotation-name my_genes \
      --gtf URL_OR_PATH --transcript-fasta URL_OR_PATH  a custom genome
  pyensembl list
  pyensembl inspect --release 93 --json
  pyensembl delete-all-files --release 75
  pyensembl prune --dry-run
```
