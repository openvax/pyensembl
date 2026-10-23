# Rebase and finish CLI error handling (#387)

## Current relevance

Current main already uses concise help and argparse errors (exit 2) for pure
selection errors, supports transcript-only FASTA input, and safely handles
explicit deletion of uninstalled genomes. Preserve those behaviors rather
than restoring obsolete branches from the old PR.

Still relevant: missing local annotation files produce tracebacks, and explicit
empty source arguments can fall through to default Ensembl selection. Keep
real parser, database and integrity failures visible with tracebacks.

## Implementation

- Rebase the existing cli/one-line-errors branch onto the newly shipped 2.32.0
  main. Resolve conflicts in favor of current release/assembly/integrity/cache
  behavior, retaining the old PR's narrow missing-data exception concept.
- Add MissingGenomeDataFile as a ValueError subclass only at the existing
  missing-remote-data conversion. Preserve library exception compatibility.
- Handle MissingLocalFile and MissingGenomeDataFile narrowly around CLI actions,
  printing `pyensembl: error: ...` and exiting 1. Missing local paths use a
  clear `file not found` message. Do not catch generic runtime ValueError,
  filesystem errors, GTF parser failures, SQL errors, or integrity errors.
- Reject explicitly empty --gtf, --transcript-fasta, --protein-fasta and
  --genome-fasta-path values during selection, before filesystem/network work.
  Keep whitespace in legitimate file paths and retain omitted optional sources.
- Pure argument errors continue to use argparse.error (exit 2). No duplicate
  CommandLineError wrapper is needed now that selection already has a narrow
  parser-error boundary. Preserve the new available/verify/latest actions.
- Close confirmed empty-source fallback issue #494 from #387. Check
  the historical corrupt-DB review claim before filing anything new; do not
  broaden this PR into parser/index recovery.
- Bump to 2.32.1; update release notes and rewrite the existing PR title/body
  around the final relevant behavior and modern exit-code convention.

## Verification and release

Use subprocesses running the real pyensembl console entrypoint to exercise
missing GTF/FASTA/DNA, unknown releases/species, invalid combinations, empty
arguments, transcript-only installs, and the current deletion behavior. Prove
empty args never reach selection/download, and malformed GTF/DB/checksum errors
retain their diagnostics. In-process CLI tests stub logging configuration.

Run targeted tests serially as well as mandatory ./lint.sh and ./test.sh, source
archive checks and docs. Preserve existing cached files exactly. Push the rebased
branch with force-with-lease against the observed remote head, pass CI, merge
PR #387, deploy from clean main, and verify the public PyPI artifacts. Read the
remaining issue queue after shipping; do not start unrelated changes.

## Verified results

- Empty GTF and mirror values selected EnsemblRelease 116 in published 2.31.2;
  confirmed with selection only and no cache creation. Fixed in this PR (#494).
- Thirty serial CLI tests pass; full lint/test suite passes 1,056 tests.
- Documentation checks pass for 26 pages / 2,187 links and 16 examples.
- The historical corrupt-DB overwrite report remains reproducible directly in
  datacache 1.29.2. Filed openvax/datacache#132; no cache-unlink workaround is
  included here. Malformed data and integrity failures retain tracebacks.
- Initial rebase used #493's final tested branch while its CI was running.
  Before pushing #387, rebase onto #493's actual merged main and publish 2.32.0.
