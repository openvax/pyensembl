# Default to the species' main genetic code (#490)

## Correction and contract

The 2.31.0 implementation unnecessarily returns an unresolved assignment for
an unrecognized contig name, even when the species' nuclear code is known.
An unplaced scaffold is an ordinary input and should use the species' main
(nuclear) table unless an organelle assignment or annotation says otherwise.
This corrects the unknown-contig policy in `genetic-code-metadata-plan.md`.

Selection order:

1. A transcript/CDS `transl_table` or `translation_table` declaration wins.
2. An explicit requested compartment or reference contig mapping selects the
   compartment. Retain validation of conflicts with established mappings and
   recognized chromosome/organelle aliases.
3. Recognized mitochondrial/plastid aliases select their respective tables.
4. Otherwise select the nuclear table, including arbitrary scaffold names and
   `genome.genetic_code()` without a contig. This is a species-specific default,
   not an unconditional table 1 fallback. For example Candida uses table 12.

Keep compartment identification separate from the final nuclear default so
`genetic_code("scaffold_42", compartment="mitochondrial")` stays valid.
`GeneticCodeAssignment.compartment` describes the selected compartment,
including the nuclear default; it does not prove assembly localization.
Retain existing table provenance. No additional API field or schema is needed.

Unknown species/tables and explicit `translation_tables={"nuclear": None}`
still produce a reasoned unresolved result. An identified organelle with a
missing table must not fall through to the nuclear table. Existing source
declaration validation and missing-GTF behavior remain unchanged.

## Implementation and compatibility

- Apply the default in the common assignment resolver used by Genome and
  Transcript; update API docstrings, the guide, and release notes.
- Preserve metadata records, serialized formats, source URLs, reference names,
  cache keys and database schema. Read existing indexes in place, unchanged.
- Add a project lesson reflecting the user's correction.
- File the PyEnsembl regression and narrow Varcode #623: its human scaffold
  table 1 fallback is intended behavior; the only remaining report concerns
  bypassing explicitly unset reference codes. Do not change Varcode here.
- Bump to 2.31.1; use a feature branch, PR, and clean-main deploy.

## Verification

- Demonstrate the new fallback tests fail against 2.31.0 before fixing code.
- Cover omitted contigs, arbitrary names, standard and nonstandard species
  codes, requested compartments, organelle aliases/mappings, source overrides,
  unresolved taxonomy, explicitly unset codes, JSON/pickle/mirrors, and both
  existing-index and GTF-stream transcript lookup.
- Run `./lint.sh` and the complete `./test.sh`, including actual read-only old
  cache fixtures and the release/cache argument matrix, in temporary storage.
- Check docs and a published Varcode consumer smoke test in the isolated env.
- Merge after required checks, run `./deploy.sh` from clean main, and verify
  the public wheel's behavior. Preserve the user's checkout, shared Python
  environment, and installed caches throughout.

## Verification results

- Before the resolver correction, the targeted regressions produced 23 failures
  against 2.31.0; after it, all 93 genetic-code tests passed.
- `./lint.sh` and `./test.sh` pass: 991 tests, 93% coverage, including installed
  read-only legacy caches, source overrides and the release/cache argument matrix.
- Published Varcode 10.18.2 selects the same scaffold and mitochondrial codes
  for human, dog, yeast, Candida and fly. Its remaining explicit-None override
  case was reproduced independently and retained as the narrower #623 report.
- Strict documentation build, 26-page/2,159-link validation, 16 existing
  executable examples and both offline genetic-code guide examples pass.
