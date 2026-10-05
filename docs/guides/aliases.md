# Look up gene name aliases

Exact-name lookup uses the symbols in your annotation. To search aliases and
previous human symbols too, load the [HGNC complete set](https://www.genenames.org/download/)
with one call. PyEnsembl downloads it once and reuses the cached file:

```python
from pyensembl import EnsemblRelease, GeneNameAliases

data = EnsemblRelease(93, species="human")  # install this release first
aliases = GeneNameAliases.download_hgnc()
genes = data.genes_by_name("p53", aliases=aliases)
print([gene.name for gene in genes])
data.close()
```

```text
['TP53']
```

Names are case-sensitive, including HGNC's lowercase `p53` alias. Alias data
is human-only and comes from HGNC independently of Ensembl. Only genes present in
the annotation are returned, with their original IDs. Ambiguous aliases return
all matching genes, including exact-name matches; do not assume the first one
is the intended locus.
`gene_ids_of_gene_name(..., aliases=aliases)` returns just the IDs.

## Reuse or select a snapshot

Cached aliases work offline. They are stored under
`aliases/homo_sapiens/hgnc` in the [PyEnsembl cache](cache.md#cache-location).
Use `cache_directory_path=` to choose another directory, or
`overwrite=True` to refresh the current complete set explicitly.

For reproducible work, keep the snapshot with your analysis. Current HGNC names
can differ from an older Ensembl release. Load a selected local snapshot with
`GeneNameAliases.from_hgnc("hgnc_complete_set.txt")`, or pass an
[archived HGNC TSV URL](https://www.genenames.org/download/)
as `source_url=` to `download_hgnc()`. Plain and gzip TSV files are supported.

## Other species or custom names

Pass your own name-to-gene-ID mapping, for example
`data.genes_by_name("old_name", aliases={"old_name": ["gene_id"]})`.
