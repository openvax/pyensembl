# Finish the unified release view (#447)

## Already shipped

Numbered and dated EnsemblRelease construction, inferred/explicit assembly
identity, dated discovery, API UUID/accession selection, old-assembly dated
selection, installed dated listings, and shared dated DNA are implemented.
EnsemblAnnotation has been removed. Retain all those APIs and cache contracts.
The remaining work concerns ordering/evidence, selection, and integrity.

## One release catalogue and explicit latest policy

- Expose one metadata view of numbered and dated releases per species/assembly,
  reusing the existing assembly catalogue and dated listing cache. It must show
  release, accession, provider, coverage, relationship evidence and selection
  status without introducing another genome/annotation constructor.
- Compare actual upstream GTFs for human numbered/dated pairs (110, 114 subset,
  115, 116 and GRCh37/75), mouse, dog, yeast, fly and Arabidopsis. Capture source
  digests, headers, gene/transcript ID.version sets and order-independent GTF
  fingerprints. Label matching identifiers separately from matching annotation
  contents. Prove subset relationships rather than inferring them from counts.
- Keep provenance and a reproducible capture/comparison tool. Tiny independent
  source excerpts exercise parsing offline; full-capture summaries describe
  precisely the evidence obtained from the complete source files.
- Add `EnsemblRelease("latest", ...)` and CLI `--release latest` as opt-in
  selection. Resolve to the latest verified selection for the requested scope,
  prefer frozen numbered data for verified ties, and serialize only the concrete
  resolved release plus selection evidence. Ordinary construction stays offline.
- Never order a date against an integer or treat every dated dataset as newer.
  Newly discovered, unranked datasets are visible and explicitly selectable,
  but cannot silently replace a verified selection. Explain the distinction
  between latest verified and the current platform dataset selected explicitly
  through the existing new-API discovery method.
- Unknown/custom assemblies can select their known dated timeline; ambiguous
  providers/assemblies need an explicit choice. Species defaults keep their
  current assembly, so tomato SL4.0 cannot roll back to dated SL3.0. Scope
  comparisons by accession, provider and annotation coverage.
- Filtered `pyensembl available` shows the unified view; keep the existing
  all-species overview. Refresh is explicit, and cached/bundled information
  works offline. Add JSON output for this metadata view if it improves the
  shared Python/CLI implementation.
- Let reference-name lookup consider installed dated releases without
  downloading or rewriting anything. Preserve preference for existing numbered
  installations unless an explicit latest policy selects otherwise. A dated-only
  installation should be reusable, including older assemblies and known explicit
  accession/provider cache layouts.

## Integrity and annotation coverage

- For new canonical dated annotation downloads, fetch the provider's md5sum.txt
  and verify each GTF/cDNA/peptide before publication. Persist the source MD5,
  manifest identity and provenance for newly acquired sources. Use datacache's
  atomic validation/publication; preserve existing source-key paths.
- Existing complete caches take their old offline, read-only fast path: no
  manifest requests, metadata migration, sidecar writes, reindex or reinstall.
  Existing mirror behavior stays valid. Legacy sources lacking acquisition
  digests must be described as unrecorded, not corrupt or silently rewritten.
- Add an explicit read-only verification API/CLI for local file integrity and
  optional comparison with current upstream manifests. Distinguish local
  corruption, missing files and changed upstream checksums. No automatic repair.
  Serialize pinned checksums so a later download cannot silently change a pinned
  dataset. Caller-requested overwrite remains explicit.
- Warn clearly for known subset annotations and expose that fact in discovery
  metadata. Keep source declarations and source coverage authoritative.
- Add a real GTF regression with coordinate-sorted/interleaved feature rows,
  including exon/CDS rows before transcript/gene rows, verifying query results
  against grouped input. Do not write a test that only checks configuration.

## Compatibility and verification

- Keep all existing reference labels, paths, files, indexes and schema versions.
  Do not alias or merge numbered and dated caches even when contents agree.
- Exercise numbered/dated selection across species, older assemblies, patched
  names and accessions; scoped latest choices; unranked/subset/mismatched
  coverage; updated listing metadata; API records; offline/read-only caches;
  JSON/pickle pinning; CLI argument combinations; and installed dated lookup.
- Use a local HTTP fixture for checksum success, missing/ambiguous manifests,
  corrupt transfers, upstream revisions and retry/concurrency behavior. Prove
  failures preserve already installed files and indexes byte-for-byte.
- Capture real source comparisons and install/query tiny data through the new
  selection and integrity routes. Re-run actual old-cache fixtures and the
  existing argument/cache matrix.
- Add any new fixtures to MANIFEST.in, run lint/full tests, source-distribution
  checks, strict docs/link/example checks, and Python 3.11–3.14 CI.
- Reconcile #447's stale description with the implementation and link existing
  downstream tracking issues; file any newly confirmed problems separately.
  Optional website-release-label aliases remain outside this scope.
- Bump to 2.32.0, merge the elaborated PR after validation, deploy from clean
  main, verify public artifacts, then rebase and finish PR #387.

## Evidence and verification results

- Complete-source comparison: human 110/115/116, mouse 116, and yeast 116 have
  identical normalized GTF records to the scoped dated counterparts. Dog 116
  and Arabidopsis 63 match identifier sets but differ in GTF records. Human 114
  dated data omit 158,710 transcript IDs; fly differs. GRCh37/75 also differs,
  and missing numbered ID versions are explicitly distinguished from bare IDs.
- All 60 dated datasets in the 52-record bundled catalogue publish checksums
  for the selected GTF, cDNA and peptide filenames (live audit 2026-10-23).
- Normal latest selection prefers the newest frozen numbered release. Explicit
  provider/coverage options restrict selection to the dated platform; within
  one accession/provider the newest known non-subset date is selected. A new
  uncompared listing cannot displace a frozen numbered release.
- 1,034 full-suite tests pass. The source archive independently checks all 133
  tracked test inputs and passes 269 offline tests, including real localhost
  transfers. Strict docs check 26 pages/2,184 links; 16 tutorial examples pass.
- Integrity tests cover concurrent acquisition, corrupt transfers, changed
  upstream/local bytes, explicit refresh, serialized pins, and exact legacy
  cache preservation. Sorted-GTF tests compare real multi-isoform queries.
- No cache schema, source path, reference label, or existing default changes.
  Numbered and mirror acquisition retain existing behavior. MD5 verification
  is for accidental corruption, and never claims sequence equivalence from GTFs.
- Varcode's dated reference-string issue openvax/varcode#616 is closed.
  Vaxrank dated CTA provenance remains tracked in openvax/vaxrank#662; it is
  downstream work, not a reason to keep two PyEnsembl constructors.
