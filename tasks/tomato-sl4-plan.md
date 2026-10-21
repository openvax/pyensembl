# Tomato SL4.0 without a separate ncRNA FASTA (#453)

## Contract

Add numbered Ensembl Genomes 63 tomato selection using the renamed species
directory and SL4.0. Declare absent ncRNA files in species/assembly metadata;
never turn arbitrary HTTP 404 errors into a successful partial installation.
Keep all previously supported source URLs, labels, serialization and cache
paths unchanged. In particular, the registered dated tomato dataset still
uses SL3.0; adding SL4.0 to numbered releases cannot relabel that dataset.

## Implementation

1. Verify upstream GTF, cDNA, peptide and DNA listings and capture a small real
   tomato annotation/sequence fixture, including independent header metadata.
2. Add an optional assemblies_without_ncrna setting to Species and consult it
   when assembling numbered transcript sources. Include it in Species equality
   and hashing so custom source selections cannot collide in cached().
3. Add a separately pinned dated_reference_name setting, defaulting to the
   latest numbered label for unchanged species. Use it for dated selections,
   installed-cache discovery, available output and catalogue generation. Pin
   tomato's existing dated label to SL3.0 without changing its accession/provider.
4. Register SL4.0 at release 63 and its renamed directory; append its verified
   assembly metadata without regenerating or changing prior catalogue entries.
5. Test real tiny installation/index/query with no ncRNA request; aliases,
   assembly assertions, CLI selection, mirror URLs, JSON/pickle and cached()
   identity; retain failure on a genuinely missing required file. Compare all
   existing selections with a pre-change snapshot and reopen old installed
   tomato fixtures offline/read-only with unchanged files and timestamps.
6. Verify the complete upstream GTF/CDNA transcript match and real installation
   in an isolated temporary cache. Existing user installations are untouched.

## Delivery

Document the declared optional source and dated/numbered distinction, bump to
2.30.1, run lint.sh and test.sh, strict docs and package checks, and open a PR
closing #453. After CI succeeds on Python 3.11–3.14, merge and run deploy.sh
from clean main. Verify the public artifacts and release tag.

## Verification results

- Full upstream installation: 33,562 genes/transcripts; every GTF transcript
  ID exactly matches a cDNA header, and every transcript sequence resolves.
- Compared 1,359 pre-change numbered/dated selections: source URLs and labels
  are unchanged. Captured old 2.30.0 tomato caches reopen offline/read-only.
- lint.sh and all 898 tests pass; strict documentation build/link checks pass.

## Separate follow-up (#481)

Build offline, provenance-bearing genetic-code metadata from captured NCBI
taxonomy. Source transcript/CDS translation tables take precedence; validated
species/strain taxonomy supplies nuclear, mitochondrial and plastid defaults.
Explicit reference compartment metadata overrides documented contig aliases;
unknown compartments and assignments remain unresolved. Preserve existing
annotation databases and serialized objects. Write a detailed API/data/test
specification after inspecting the source attributes and downstream contract,
before changing that implementation. Deliver it in its own minor-version PR.
