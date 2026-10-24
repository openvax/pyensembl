# PyEnsembl #456: remove typechecks and establish focused static checking

## Compatibility contract

Replace the three database typechecks calls and chromosome predicates with plain
isinstance checks. Keep invalid-type TypeError and empty-column ValueError.
Builtin int subclasses (including bool) remain accepted and numpy integers stay
unsupported. Keep all accepted chromosome normalization, locus SQL query behavior,
return values, defaults, optional coordinates, source metadata, file paths and
cache reuse. Do not migrate, reinstall, rename, write or reindex user caches.

A separately filed bug discovered during baseline review allows invalid numeric
chromosomes through a warm memoization cache. Validate type before lookup, making
rejection independent of call history; retain memoization for accepted inputs.

## Implementation

Remove the runtime dependency and imports. Annotate normalize_chromosome,
normalize_strand and the guarded locus query methods plus a small directly
supporting surface. SQLite result cells are dynamic, so return annotations must
not falsely claim every queried column is text. Preserve arbitrary column values.
Keep the typing effort bounded: start mypy with the changed modules, annotate
state required for real checking, and skip unrelated imported modules rather
than blanket-disable errors. Don't add another runtime validation package.

Ship py.typed with wheel/sdist and document the partial annotation scope. Add
mypy as a development dependency, invoked by lint.sh, which CI already runs.
Test typing with a real small consumer and reveal/assert types; prove an invalid
consumer is rejected. Add runtime regressions against a tiny temporary GTF/DB
for str/int subclasses, bool, invalid types, None, floats, numpy values, empty
columns, and successful string/numeric-column locus queries. Include cold/warm
normalization behavior and unchanged installed-cache snapshots.

## Verification and release

Work from fresh main after the earlier authorized releases, bump the next patch
version, update release notes and add new offline tests to check_sdist.py. Run
lint.sh (ruff+mypy), test.sh, focused serial checks, extracted archive tests,
artifact/marker checks, docs, pip check and Python 3.11–3.14 CI. Merge through PR,
run deploy.sh from clean main and verify the public artifacts. Shared environment
and the user's modified pyproject.toml/tasks/lessons.md are left intact.
