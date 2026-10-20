# Dataset constructors

`EnsemblRelease` selects numbered and dated Ensembl datasets, infers assemblies,
and supports explicit assembly/provider choices. See
[assembly selection](../guides/assembly-selection.md) for examples, metadata,
API discovery, and cache compatibility. Both it and custom `Genome` objects
provide the shared [Genome methods](genome.md).

::: pyensembl.EnsemblRelease
    options:
      inherited_members: false

## Compatibility constructor

`EnsemblAnnotation` retains its historical signature, imports, serialization,
and cache defaults. New code should use `EnsemblRelease`.

::: pyensembl.EnsemblAnnotation
    options:
      inherited_members: false
