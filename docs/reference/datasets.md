# Dataset constructors

`EnsemblRelease` selects a numbered [Ensembl](https://ftp.ensembl.org/pub/)
release, such as 116, or a dated release of [the new platform](../ensembl-platform.md),
such as `"2026_04"`. `EnsemblAnnotation` selects any assembly accession,
provider and dated geneset from the new platform. Both provide the shared
[Genome methods](genome.md).

::: pyensembl.EnsemblRelease
    options:
      inherited_members: false

::: pyensembl.EnsemblAnnotation
    options:
      inherited_members: false
