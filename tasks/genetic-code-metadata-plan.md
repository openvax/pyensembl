# Offline source and taxonomy genetic-code metadata (#481)

Policy correction in 2.31.1: the arbitrary-contig unknown behavior below is
superseded by [the species-default plan](species-genetic-code-default-plan.md).
Unplaced scaffolds and no-contig lookups now use the species' nuclear code.
The original Varcode #623 scaffold example was intended fallback behavior;
that issue now covers only explicitly unset reference tables.

## Scope and invariants

Expose translation-table assignments and provenance from PyEnsembl, without
implementing translation, modifying sequence data, or guessing a code from
an organelle name alone. Preserve old cache paths, source URLs, database
schema compatibility, serialized objects and installed indexes. All runtime
lookups are offline. Use isolated development/release environments.

## Public contract

- `Genome.reference_metadata`: a copy of a serializable pinned metadata
  record. One optional constructor argument with that name allows custom
  references to supply taxon identity, compartment assignments, and explicit
  translation-table overrides. EnsemblRelease infers the record normally.
- `Genome.taxon_id`, `species_taxon_id`, `nuclear_translation_table`,
  `mitochondrial_translation_table`, `plastid_translation_table`, and
  `contig_compartments` expose the resolved reference metadata.
- `Genome.genetic_code(contig=None, *, compartment=None)` returns one immutable
  provenance record: table_id, compartment, source, taxon_id, source_url and
  an explicit reason when unresolved. A compartment can be requested without
  a contig. Supplied inconsistent constraints raise instead of guessing.
- `Transcript.source_translation_table` exposes only annotation declarations;
  `Transcript.genetic_code` resolves the source declaration before the
  reference defaults; `Transcript.translation_table` exposes the resulting
  table ID or None. Conflicting declarations on transcript/CDS rows raise.

## Taxonomy and provenance

Capture official Ensembl species/taxonomy mappings and NCBI Taxonomy XML with
lineage/parent records, synonyms, nuclear/mitochondrial IDs, and plastid pgcode.
Capture NCBI gc.prt solely to validate table IDs/names; do not duplicate a
translation engine. Generate and package a small reviewed JSON snapshot,
with retrieval date, source URLs and content digests. Generation consumes
captured files and is deterministic; capture is a separate explicit command.

Resolve registered species and exact assembly selections, including dog
breeds and mouse strains, from verified mappings. Distinguish source taxon
from its species ancestor. Accept a parent/child relationship only with
captured lineage evidence; unrelated supplied identities are errors. Unknown
taxa and unavailable code assignments remain explicit. Retain actual and
species taxon IDs from new Ensembl API responses and validate their agreement.

## Compartment resolution

Explicit per-contig reference metadata takes precedence over documented
aliases. Recognize mitochondrial/plastid aliases and a conservative set of
nuclear chromosome names; an arbitrary remaining contig is unknown, not
automatically nuclear. A taxon/compartment assignment still must exist before
returning a default code. Do not infer taxonomy from filename fragments or
custom display labels. Exact registered reference labels retain legacy use.

## Annotation and cache compatibility

New annotation databases retain optional transl_table/translation_table
columns, including an empty value when absent, without changing the schema
version or rebuilding existing databases. For old indexes lacking those
columns, lazily inspect the already installed local GTF and cache its sparse
declarations only in memory. No network, database migration, new sidecar, or
reindex on lookup. A missing required local annotation source must be reported,
not treated as proof that there was no source override. Source fingerprint
changes invalidate any in-memory metadata. Preserve ordinary index/query
behavior and keep the new lookup separate from download-triggering db access.

Review refinement: new indexes also record the GTF stat fingerprint for these
optional attributes. Use indexed declarations only when that fingerprint
still matches and both transcript/CDS features were examined. Otherwise read
the installed GTF. This prevents a modified remote-source file or a restricted
custom index from falsely reporting no override. Existing indexes are never
changed to add the fingerprint. Accept the same quoted/repeated attribute
forms as the main GTF parser, including quoted semicolons in unrelated fields.

Pin metadata through Genome/EnsemblRelease JSON, pickle and mirrors. Old
saved objects gain inferred metadata in memory only. Metadata must not alter
source or cache identities; supplied semantic overrides must remain distinct
in in-memory object caching.

## Acceptance and release

Test actual GTF/CDS overrides, absent attributes, invalid IDs and conflicting
rows; human/mouse mitochondrial 2, yeast 3, fly/worm 5, Candida nuclear 12,
plant plastid 11, parent/strain resolution and conflicting taxa. Test unknown
contigs/taxa, explicit compartment overrides, API taxon propagation, mirrors,
serialization and snapshot preservation across catalogue changes. Reopen
old cached databases read-only and prove no downloads, writes or rebuilds.
Use source-derived examples independently of generated config, and verify
the snapshot generator against captured NCBI metadata.

Run lint.sh/test.sh, strict docs/tutorial and installed-wheel tests, plus a
focused read-only consumer smoke test against Varcode. File any consumer
integration defects; this task changes only PyEnsembl. Ship a separate 2.31.0
PR after #453, pass Python 3.11–3.14 CI, merge, run deploy.sh from clean main,
and verify public artifacts and the release tag.

## Implementation verification

- `./lint.sh` passed; `./test.sh` passed 983 tests, with 93% coverage.
- Strict documentation build, 26-page/2,159-link audit, and all 16 executable
  documentation examples passed. Wheel and sdist build successfully.
- The installed wheel passed 167 metadata, source-declaration, tomato and
  release/cache combination tests, imported from an isolated target directory.
- An actual 2.30.1 index fixture recovers a formerly unindexed source override;
  actual 2.30.0 tomato indexes recover absent declarations and taxonomy defaults.
  Read-only before/after snapshots preserve all files, bytes, modes, inodes and
  timestamps. The existing release/cache argument matrix also passes.
- Independent primary captures reproduce all 43 mappings and 42 NCBI records.
- Published Varcode 10.18.2 consumes the expected human/mouse, yeast, fly/worm,
  Candida and plant assignments and a GTF override. It still falls back to
  nuclear table 1 after an explicit unknown compartment; filed downstream
  [Varcode #623](https://github.com/openvax/varcode/issues/623). No Varcode or
  shared-environment files were changed.
- All work, caches, builds and environments are isolated under `/private/tmp`.
