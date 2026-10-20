# Remove the unused constructor and exercise installed cache selections

## Intent and authorization

The user explicitly approved removing EnsemblAnnotation in the next version,
superseding #482's indefinite compatibility-shim requirement. They also asked
for real locally cached data created with different EnsemblRelease argument
combinations. Keep the existing file-cache preservation requirement unchanged.

Ship this as 2.30.0 through a feature PR, then merge and run deploy.sh from clean
main. Keep the original user checkout and shared Python environment untouched.

## Product changes

- Delete the EnsemblAnnotation implementation/module, package export, and
  constructor-only branches in EnsemblRelease. Do not add an alias, hidden
  compatibility class, warning wrapper, or serialized-object adapter.
- Keep EnsemblRelease's numbered and dated behavior, explicit accession cache
  layout, reference names, URL spelling, provider/coverage selection, and
  metadata serialization unchanged.
- Keep old EnsemblRelease JSON/pickle loading. Old EnsemblAnnotation imports and
  serialized class paths intentionally cease to work; explain the replacement
  in one short removal note. Existing cache files remain reusable by an
  equivalent EnsemblRelease call, with no rewriting or regeneration.
- Update live documentation and move the actual assembly/date integration tests
  to the unified constructor. Preserve historical fixtures as independent
  evidence instead of regenerating them with the implementation under test.

## Installed-data matrix

Use a new isolated temporary cache and tiny real GTF/FASTA data. Exercise actual
acquisition, decompression, database creation, sequence indexing, and queries;
assert resulting file paths and contents, not just constructor attributes.
Transport fixtures must not contact Ensembl or affect users' installed files.

Cover:

1. Numbered human release 75 via default species, common/Latin names, build,
   exact patch name, and accession assertions: same installed annotation.
2. Dated human default via inferred species, broad/exact names, and explicit
   default coverage: same species-layout data including alternate loci.
3. Explicit accession via either keyword, normalized spelling, optional species
   and provider: same accession-layout primary data; preserve the intentional
   distinction from the species default.
4. Historical dated GRCh37 inference versus build/name/accession selection:
   identical selected data in the explicit layout, with the existing mirror
   and previously installed-selection preservation rules respected.
5. Dog reference and breed selections: synonyms and exact selectors reuse the
   intended cache; dog/breed assemblies and providers do not cross-contaminate.
6. Different dates, providers, coverage, and assembly accessions: cache choices
   remain separate even when source basenames match. Use distinct synthetic
   gene/transcript/protein/DNA contents to catch wrong-cache reuse.
7. Custom labels/directories, remote versus local DNA, and masked DNA: reopen
   the same installed data through constructor, JSON, pickle and install CLI
   arguments where applicable. Keep DNA attachment independent of annotation.
8. Invalid/mismatched species, patch/accession, and ambiguous breed dates: fail
   without writes. Repeat installations and read-only/offline reopening must
   not download, parse, rewrite, or change existing cache bytes, permissions,
   inodes or modification/change timestamps.
9. Published 2.28.3 fixture caches: retain old Release deserialization and reopen
   old Annotation-created data using the equivalent EnsemblRelease arguments.

Keep the matrix understandable by grouping equivalent selections and checking
cross-group differences. Extend existing DNA-sharing tests rather than
reimplementing the shared DNA cache. Keep fixtures small and deterministic.

## Audit and validation

Search live downstream source/tests for removed imports, excluding vendored
copies of PyEnsembl itself. File concrete discovered defects in their owning
repository and link them from the PR; fix cache defects within this scope.

Run focused installed-data tests, serial CLI tests, lint.sh, test.sh, strict docs
and tutorial examples. Build wheel/sdist, verify the removed module is absent,
and run the new matrix against the installed wheel. Check final CI on Python
3.11–3.14, merge, run deploy.sh, verify public artifacts and release tag, then
review the next relevant dependency block without starting unrelated changes.

## Matrix finding and revised implementation

The local/remote/masked DNA matrix reproduced mutation of genome_fasta.json
when an older indexed local source is reopened. Recognize existing private
source artifacts or an existing matching manifest without writing migration
metadata. Record source-specific evidence only for a new download-only local
attachment with no owned cached artifacts, so subsequent reopening can avoid
changing the advisory manifest after other attachments are installed. Expand
regressions to two local inputs, compressed local DNA, download-only attachment
recording, and read-only reuse, then rerun the full validation gates.

## Local verification

- Filed removal/matrix issue #484 and the discovered manifest defect as #485.
- All 67 constructor combinations across 17 selections create fresh caches,
  then reopen through equivalent arguments, cached(), JSON, pickle and CLI
  selection without changes. Additional tests cover coexistence, six optional
  DNA attachments, invalid constraints, and API accession/UUID discovery for
  human and dog.
- Historical 2.28.3 fixture files remain unchanged and reopen through the
  supported constructor; old EnsemblRelease saved objects still load.
- lint.sh passes; test.sh passes all 889 tests (93% coverage). Serial CLI/matrix
  checks passed, as did 195 installed-wheel checks, strict documentation/link
  validation, 16 real-data tutorial examples, and pip check.
- Both build artifacts include the assembly catalogue and omit the removed
  module. Downstream search found no live imports in the relevant consumers.
- CI, merge, deployment, and public-artifact verification remain release gates.
