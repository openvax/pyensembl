# Self-contained source-distribution tests (#489)

The published sdist includes test modules but omits their supporting Python
modules and all 48 tracked fixtures. Ship the required test package, fixtures,
and lint/test entry points so an extracted source archive can run its suite.
Ordinary wheel contents and runtime/cache behavior remain unchanged.

The archive also needs the root contributor/release documents referenced by
its included documentation pages, including the new artifact-test instructions.

Use explicit fixture paths in MANIFEST.in: recursive inclusion would pick up
untracked developer data and generated indexes. Keep archived legacy-cache
fixtures, whose gzip archives are intentional test inputs, while excluding
loose generated databases, pickles and indexes. New fixtures must be added to
the manifest; verification compares the resulting archive with tracked source
fixtures independently of the manifest's implementation.

Add an artifact verifier that extracts a built sdist into a temporary directory
and executes offline tests covering FASTA/GTF parsing, assembly metadata,
genetic-code captures and old cache fixtures. Assert imports come from the
extracted package and use temporary cache roots. When run from a Git checkout,
also compare tracked fixture names and bytes with the archive, rejecting extra
fixture files. It must remain usable without Git in an extracted sdist.

Run the verifier in a separate CI job and before upload in deploy.sh. Document
the command. First demonstrate failure on the existing public 2.31.1 archive;
then build with sentinel generated/local files present, verify their exclusion,
and run the full lint/test suite from both checkout and extracted archive.
Bump to 2.31.2, merge through a PR, deploy from clean main and verify public
artifacts. Preserve the user's checkout, installed caches and shared environment.

## Verification

- The verifier rejected public 2.31.1: 48 missing fixtures and three missing
  helper modules. The corrected archive matches all 125 tracked test inputs.
- Six untracked sentinel files (GTF, JSON, SQLite, pickle, FAI and GZI) were
  excluded from the archive despite being present at build time.
- Extracted-archive offline checks: 226 passed, importing the extracted package.
- `./lint.sh` and `./test.sh`: 991 passed in the checkout and independently in
  the extracted source archive, using isolated copies of required genome data.
- Strict documentation build and 2,159-link audit passed.
