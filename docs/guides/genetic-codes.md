# Find a reference or transcript's genetic code

PyEnsembl exposes NCBI translation-table IDs and their provenance. It does
not translate sequences. Use a transcript's assignment when available: a
declaration in its annotation takes precedence over the reference default.

## Reference defaults without downloads

```python
from pyensembl import EnsemblRelease

yeast = EnsemblRelease(116, species="yeast")
assert yeast.taxon_id == 559292           # the S288C reference strain
assert yeast.species_taxon_id == 4932     # Saccharomyces cerevisiae
assert yeast.genetic_code("MT").table_id == 3
assert yeast.genetic_code(compartment="nuclear").table_id == 1

plant = EnsemblRelease(63, species="tomato")
assert plant.plastid_translation_table == 11
```

`nuclear_translation_table`, `mitochondrial_translation_table` and
`plastid_translation_table` return an integer or `None`. The assignments
come from a [bundled, verified taxonomy snapshot](../dev/genetic-code-catalogue.md),
with independent values for each compartment. Mitochondrial assignments
include human/mouse 2, budding yeast 3 and fly/nematode 5; Candida albicans
has nuclear table 12. A mitochondrial contig name alone does not select a
table. The reference's taxon or an explicit assignment is also needed.

`genetic_code()` returns an immutable `GeneticCodeAssignment` with `table_id`,
`compartment`, `source`, `taxon_id`, `source_url`, and `reason`. Unknown
assignments have `table_id=None` and an explanation. There is no automatic
standard-code fallback for arbitrary unplaced contigs or unknown taxa.

## Transcript source declarations

After the annotation has been [installed](cache.md), read:

```python
transcript = yeast.transcript_by_id(yeast.transcript_ids()[0])
source_only = transcript.source_translation_table
assignment = transcript.genetic_code
effective_id = transcript.translation_table
```

`source_translation_table` reads `transl_table` or `translation_table` on the
transcript and its CDS rows. It returns `None` when no declaration is present.
Repeated agreeing values are accepted; different values or invalid NCBI IDs
raise `ValueError`. A valid source declaration wins over a different taxonomy
default or explicit reference default, even on an unknown compartment.
`assignment.source` is `annotation`, `reference`, `taxonomy`, or `unresolved`.
For an annotation declaration, `source_url` identifies the configured GTF
path or URL. The compartment can still be `None` if it cannot be established.

New annotation indexes retain these optional attributes. Old indexes remain
usable unchanged: the first lookup scans the installed GTF for the sparse
declarations and keeps them only in memory. There is no download, migration,
sidecar, or reindex. If a configured GTF is missing, transcript lookup raises
instead of treating that absence as evidence that no source override exists.
Reference-only `genome.genetic_code()` does not need a GTF.

## Custom references and contigs

One optional `reference_metadata` argument configures the reference:

```python
from pyensembl import Genome

genome = Genome(
    "study-reference", "study-annotation",
    reference_metadata={
        "taxon_id": 7227,
        "contig_compartments": {"scaffold_42": "mitochondrial"},
        "source_url": "https://example.org/study-reference",
    },
)
assert genome.genetic_code("scaffold_42").table_id == 5
assert genome.genetic_code("unplaced_1").table_id is None
```

An optional `translation_tables` mapping can explicitly set `nuclear`,
`mitochondrial`, or `plastid` to a valid integer table ID. `None` explicitly
leaves that compartment unresolved. Unspecified compartments retain their
validated taxonomy defaults. Table IDs are metadata: a translation engine
may reject codes whose stop/sense behavior needs additional context.

Explicit `contig_compartments` take precedence. Matching is exact first,
then case-insensitive; conflicting case aliases raise. Fallback names are:

| Compartment | Case-insensitive aliases |
| --- | --- |
| Mitochondrial | `MT`, `M`, `chrMT`, `chrM`, `mito`, `mtDNA`, `mitochondrion`, `mitochondrial`, `mitochondrion_genome`, `mitochondrial_genome` |
| Plastid | `Pt`, `CP`, `chrC`, `chrPt`, `chloroplast`, `plastid`, `apicoplast` |
| Nuclear | Positive chromosome numbers, `X`, `Y`, `W`, `Z`, Roman numerals I–XVI, and chromosome arms 1L–4R; optionally prefixed by `chr` |

These names are conventions, not assembly evidence. Supply reference
metadata for nonstandard names. You can also call
`genome.genetic_code(compartment="plastid")` without a contig. When both are
given, an inconsistent known compartment raises rather than ignoring the
constraint.

## Taxonomy and saved objects

Registered species infer taxon identity offline. A custom `Genome` can infer
it from an exact registered `reference_name`; arbitrary labels and filename
fragments are not matched. An `EnsemblRelease` uses its selected species or
source API record, independent of its display label. `from_ensembl()` retains
the API's taxon and species-ancestor IDs and validates known conflicts.
Strain/parent differences are accepted only with captured lineage evidence.
Unknown positive taxon IDs remain unresolved; they do not trigger requests.

`genome.reference_metadata` returns a copy of the pinned record, including
source URLs, retrieval time and taxonomy digest. JSON, pickle and
`with_mirror()` preserve it. Loading a saved record does not refresh its
codes from a newer catalogue. Old saved objects infer metadata in memory.
Metadata overrides distinguish in-memory objects, while annotation source
URLs, reference labels and disk cache paths stay unchanged. CLI install hints
install the same files; custom metadata belongs to the Python object and
its serialized record, not the on-disk installation.
