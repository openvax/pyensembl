# PyEnsembl #470: reuse a mirror installation from Python

## Problem and boundary

Since 2.27.5 new downloads are keyed by their complete source URLs. The CLI
already builds a custom Genome for a flat-directory mirror; ordinary
EnsemblRelease constructs canonical URLs. They intentionally have different
cache identities. Existing basename-only installations are addressed separately
by #474; this change must preserve that compatibility contract and new-source
isolation. No copying, relabeling, or implicit registration of mirror indexes.

## API and implementation

Add `EnsemblRelease(...).with_mirror(mirror)` returning a source-explicit Genome.
Use the release's existing filename selection for GTF, transcript and protein
FASTA, and reference DNA when requested. Retain explicit local DNA; support
the CLI's existing explicit DNA override. Route `--custom-mirror` through this
same method and remove the duplicated URL rewriting from the CLI. The returned
Genome uses its existing serialization and identity contracts.

Document a complete install/reopen example, including passing the returned
Genome explicitly to Varcode. Explain that default EnsemblRelease/Varcode release
numbers require canonical installation for new caches. Mirror-aware CLI queries
must repeat the same mirror. Existing mirror downloads do not need migration:
construct the same source URLs. Include these distinctions in CLI help and
cache documentation. Do not change default acquisition or existing files.

## Verification

Use a local HTTP server serving tiny real GTF/FASTA data under independently
specified Ensembl filenames. Run the actual install CLI in a subprocess with
a fresh cache. Stop the server, then reopen/query through the Python API in
another process. Check canonical and a second mirror do not claim this install;
snapshot installed files to rule out rebuilding/copying. Exercise pickling/JSON,
mirror CLI selection, local/remote reference-DNA overrides and numbered/dated
file selection. Run existing cache isolation and CLI regressions, lint, and the
full suite. Bump the version after refreshing main to avoid conflicting with
the active #474 release. Ship only after green CI, merge, clean-main deployment
and live PyPI verification; refresh the shared installed package normally.
