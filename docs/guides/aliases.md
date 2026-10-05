# Look up gene name aliases

Exact-name lookup uses the symbols in the selected annotation. To also search
aliases and previous human symbols, download an [HGNC complete-set TSV
snapshot](https://www.genenames.org/download/) and load it explicitly:

```python
from pyensembl import EnsemblRelease, GeneNameAliases

data = EnsemblRelease(93, species="human")  # install this release first
aliases = GeneNameAliases.from_hgnc("hgnc_complete_set.txt")
genes = data.genes_by_name("p53", aliases=aliases)
```

Names are case-sensitive, including HGNC's lowercase `p53` alias. Alias data
is read locally. Keep the snapshot with your analysis: current HGNC
nomenclature can differ from an older Ensembl release. Only genes present in
the annotation are returned, with their original IDs. Ambiguous aliases return
all matching genes, including exact-name matches; do not assume the first one
is the intended locus. HGNC data is human-only. For other species or custom
annotations, pass a mapping such as `aliases={"old_name": ["gene_id"]}`.
`gene_ids_of_gene_name(..., aliases=aliases)` returns just the IDs.
