# Adopt published datacache recovery in PyEnsembl

## Scope and compatibility

Require datacache 1.29.3, the first published version containing the recovery
fix in openvax/datacache#133. An existing installation must receive that version
when PyEnsembl is upgraded. PyEnsembl must delegate recovery to datacache rather
than delete a user's SQLite cache or invent a second recovery mechanism.

Normal opens and installs retain their current behavior. Only explicit
`--overwrite` may recover a corrupt annotation database. Preserve healthy cache
reuse, source files, indexes, paths, permissions, links, ownership and rollback
semantics. Use disposable data in tests, never corrupt an installed user cache.

## Implementation and evidence

- Add a CLI integration regression using an actual tiny GTF and database.
  Corrupt the temporary database, verify normal installation reports the
  existing error without rewriting it, then invoke `install --overwrite` and
  query the reconstructed gene. Confirm the database inode and permissions
  survive and annotation source bytes are unchanged.
- Exercise recovery failure with malformed source data; preserve the old
  database and source when rebuilding cannot complete. Existing healthy-cache
  reopen snapshots continue to prove no new writes during normal reuse.
- Demonstrate the recovery regression against datacache 1.29.2 in a separate
  temporary environment, and success against published 1.29.3. Do not assert a
  hardcoded copy of the dependency declaration in a test.
- Raise the dependency floor, bump PyEnsembl to 2.32.2, and document the explicit
  recovery behavior and the upstream runtime requirements without suggesting
  migration or automatic repair.
- File the integration gap as a PyEnsembl issue and link it and datacache#132 /
  #133 in the PR. Run mandatory lint/full tests, serial CLI tests, extracted
  source-archive tests, docs and package checks. All environments and test data
  remain separate from the shared virtual environment and user checkout.

## Shipping and subsequent authorized work

Wait for Python 3.11–3.14 CI, merge and run deploy.sh from a clean main clone.
Verify public artifacts and the release tag. Then implement OncoRef #609,
finish the existing Vaxrank #665 for #662, and finally PyEnsembl #456. Each
stage gets its own versioned PR, repository-required validation and release.
