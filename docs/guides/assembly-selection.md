# Choose an assembly and annotation

Use `EnsemblRelease` for both numbered Ensembl releases and dated annotations
from the new Ensembl platform. Pin the annotation release for reproducibility.
The assembly is normally inferred from the species and selected release:

```python
from pyensembl import EnsemblRelease

human = EnsemblRelease(75)                              # GRCh37.p13
current = EnsemblRelease("2026_04", species="human")    # GRCh38.p14
dog = EnsemblRelease("2020_10", species="dog")          # ROS_Cfam_1.0
historical = EnsemblRelease("2013_09", species="human")  # GRCh37.p13
```

Construction is offline and does not create cache directories. Call
`download()` and then `index()` to install a dataset; existing installations
remain usable in place. Use [custom genomes](custom-genomes.md) for your own
matched GTF and FASTA files.

## Numbered releases

PyEnsembl supports numbered releases through Ensembl 116 and Ensembl Genomes
63. Common and Latin species names are accepted. The species/release mapping
selects the historical assembly, including its actual patch revision where
recorded by Ensembl.

```sh
pyensembl available
pyensembl install --release 75 --species human
pyensembl install --release 104 --species dog
```

Human release 75 uses GRCh37.p13; human release 93 uses GRCh38.p12. Dog release
104 uses CanFam3.1, while 105 uses ROS_Cfam_1.0. An explicit assembly assertion
must agree with the release:

```python
EnsemblRelease(75, assembly="GRCh37.p13")  # accepted
EnsemblRelease(75, assembly="GRCh38")     # raises ValueError
```

The existing `--reference-name GRCh37` shortcut selects human release 75 when
no release is supplied. `genome_for_reference_name("GRCh37")` prefers an
installed release, then a downloaded release, then the newest supported
numbered release. Installed dated releases participate, including explicit
accession/provider layouts; an installed numbered release keeps precedence.
Selection reads existing caches without changing them. Deletion commands
require a concrete release and reject `latest`.

## Dated releases

A date such as `"2026_04"` identifies an annotation directory, not a website
release. A website label such as `"2026-07"` is a different identifier and is
rejected by the ordinary constructor.

```sh
pyensembl install --release 2026_04 --species human
pyensembl install --release 2013_09 --assembly GRCh37
```

Selection preserves the registered default when it has that date. When the
bundled catalogue identifies a unique historical dataset, it can infer that
assembly instead: human `2013_09` selects GRCh37.p13. An existing cache directory,
old serialized selection, or custom mirror retains its previous selection.
Unknown dates retain the existing default selection; absence from the bundled
catalogue is not proof that a newly published date does not exist.

To discover dates for the registered default assembly:

```python
from pyensembl import available_dated_releases, available_annotation_dates

available_dated_releases("human")
available_annotation_dates("GCA_000001405.14", "ensembl")
```

These explicit discovery functions fetch and cache directory listings on first
use, then support offline reuse. `refresh=True` checks again. Ordinary dataset
construction does not fetch these listings. A missing remote file produces an
`UnpublishedDateError` with available dates when discovery confirms that the
annotation date is absent.

Dates are not a universal chronology: Arabidopsis uses `2010_09`, and fly's
dated `2022_07` annotation predates its final numbered annotation. Matching
assemblies do not make two annotations equivalent. There is no automatic
conversion between numbered and dated releases or rolling default.
For example, human `2024_11` is a GENCODE subset without non-canonical lncRNA
transcripts; numbered release 114 contains the complete annotation.

## Unified catalogue and opt-in latest

```python
from pyensembl import available_releases

rows = available_releases("human")
chosen = EnsemblRelease("latest", species="human")  # concrete release 116
older = EnsemblRelease("latest", assembly="GRCh37")  # concrete release 75
```

`available_releases` combines numbered and dated entries for one assembly.
Without an assembly selector, it uses the species' latest numbered assembly;
tomato therefore stays on SL4.0. Build names, patched names, accessions, and
species constraints follow the ordinary constructor. Each entry includes
assembly identity, provider, coverage, subset status, comparison evidence, and
whether it is the `latest` choice. Numbered entries followed by dated entries
are a display grouping, **not** a cross-platform chronological ordering.

The catalogue reads bundled metadata and any cached dated listings offline.
`refresh=True` explicitly fetches listings. Dates discovered after the bundled
comparison audit appear as uncompared; they cannot displace a frozen numbered
release merely because their date looks recent.

`latest` is opt-in and conservative: use the newest frozen numbered release
for the selected assembly unless a reviewed ordering explicitly supersedes it.
Prefer numbered data on verified ties. Supplying `provider=` or `include_alt=`
restricts selection to the dated platform; within one accession/provider it
chooses the newest known date, excluding known subsets. Ambiguous dated scopes
require an accession/provider. Existing no-argument and CLI defaults stay the
same. Use `from_ensembl` below for current platform discovery through the API.

The result stores a **concrete** `.release` and resolved metadata, including the
selection policy. JSON and pickle restore that snapshot without resolving
`latest` again. Keep that serialized object or record its concrete arguments
for reproducible analyses.

```sh
pyensembl available --species human
pyensembl available --assembly GRCh37 --json
pyensembl available --species dog --refresh
pyensembl install --release latest --species dog
```

Unfiltered `pyensembl available` retains the all-species overview. Filtered
catalogue output is offline unless `--refresh` is requested.

### What the comparisons establish

The October 2026 audit compares complete GTF files, not just their headers or
row counts. Human GRCh38 uses the including-alt GTF; other rows use the ordinary
GTF selected by `EnsemblRelease`. Equality excludes headers and row order but
includes all data columns, attributes, and duplicate rows.

| Species | Numbered ↔ dated | Verified relationship |
| --- | --- | --- |
| Human GRCh38 | 110 ↔ `2023_03`, 115 ↔ `2025_12`, 116 ↔ `2026_04` | Same normalized GTF records |
| Human GRCh38 | 114 ↔ `2024_11` | Dated transcript IDs are a strict subset; 158,710 numbered transcripts absent |
| Human GRCh37 | 75 ↔ `2013_09` | Different; the numbered source also omits ID versions |
| Mouse | 116 ↔ `2026_04` | Same normalized GTF records |
| Dog | 116 ↔ `2020_10` | Same gene/transcript ID.version sets; different GTF records |
| Yeast | 116 ↔ `2018_10` | Same normalized GTF records |
| Fly | 116 ↔ `2022_07` | Different gene/transcript sets |
| Arabidopsis | 63 ↔ `2010_09` | Same gene/transcript ID.version sets; different GTF records |

These comparisons make no claim about FASTA sequence equality. Primary-only
human selections do not inherit including-alt evidence. Even identical GTF
records do not alias or merge annotation caches. Source URLs, checksums, counts,
set differences, and reproducible capture instructions are described in the
[assembly catalogue maintenance guide](../dev/assembly-catalogue.md#annotation-comparisons).

## Dated annotation integrity

New downloads from the canonical dated platform verify each compressed GTF,
cDNA, and peptide file against its provider's `geneset/md5sum.txt` **before**
publishing the file. A receipt records the source MD5, manifest URL/digest,
and retrieval time. JSON/pickle pin recorded checksums. Downloading a pinned
source with a changed upstream MD5 fails; explicit `overwrite=True` accepts a
new acquisition after validating it against the current manifest.

Existing complete caches retain their original offline path. Normal reuse
neither fetches manifests nor adds receipts, changes paths, or rebuilds indexes.
Mirrors and numbered downloads keep their existing acquisition behavior.

```python
reports = current.verify_sources()                   # read local files only
reports = current.verify_sources(check_remote=True)  # also read upstream manifests
```

```sh
pyensembl verify --release 2026_04 --json
pyensembl verify --release 2026_04 --check-remote
```

Verification never repairs or rewrites data. Each source reports `verified`,
`changed`, `unrecorded`, or an availability failure such as `missing`.
`unrecorded` means no checksum was captured when the source was acquired; it
is not a corruption finding. Optional upstream results separately distinguish
an unchanged manifest entry from a revised one. The command exits with status
1 for a missing/unreadable/changed source or changed upstream checksum, and 0
for verified or unrecorded sources. MD5 detects accidental changes; it is not
cryptographic proof of authenticity. These checks cover annotation sources;
reference DNA has its separate checksum/index verification.

Known subsets such as human `2024_11`, and API annotation labels ending in
`_subset`, produce a warning. `annotation_is_subset` exposes that status.
`annotation_label` preserves the provider's API label when supplied; it is
separate from `.release` and the existing `.annotation_version` compatibility
field, which both identify the selected PyEnsembl release.

## Explicit assemblies and dog breeds <a id="choose-the-dataset-explicitly"></a>

Use an assembly argument when the requested date is ambiguous, when selecting
a nondefault assembly, or to assert the assembly your analysis requires.
All three forms below identify the same GRCh37.p13 assembly:

```python
EnsemblRelease("2013_09", assembly="GRCh37")
EnsemblRelease("2013_09", assembly="GRCh37.p13")
EnsemblRelease("2013_09", assembly="GCA_000001405.14")
```

A broad name resolves the patch appropriate to the dataset. It does not mean
the newest patch in every historical release. Conflicting species, build,
patch, and accession arguments raise `ValueError` before acquisition.

Dog breeds remain distinct assemblies. Basenji and Great Dane both publish an
annotation dated `2019_08`, so `EnsemblRelease("2019_08", species="dog")` asks
for a more specific selection:

```python
basenji = EnsemblRelease("2019_08", assembly="Basenji_breed-1.1")
great_dane = EnsemblRelease("2019_08", assembly="UMICH_Zoey_3.1")
```

Species is inferred from known exact selections. Registered breed aliases
also work, such as `species="dog_basenji"`. An unregistered accession remains
usable, but does not acquire human identity merely because human is the
ordinary constructor default.

## Assembly metadata

These read-only attributes come from one resolved metadata record, regardless
of whether the assembly was inferred or supplied explicitly:

| Attribute | Meaning | Human release 75 |
| --- | --- | --- |
| `assembly` | Build/coordinate-system identifier | `GRCh37` |
| `assembly_name` | Complete official name, including patch when present | `GRCh37.p13` |
| `assembly_accession` | Exact versioned GCA/GCF accession | `GCA_000001405.14` |
| `reference_name` | Existing display/cache compatibility label | `GRCh37` |
| `release` | Concrete annotation release | `75` |

Names such as `CanFam3.1` and `ROS_Cfam_1.0` retain their meaningful version
suffixes. `reference_name` is deliberately independent of the official name:
for example, the old fly cache label remains `BDGP6.54` while the dated
platform calls the assembly `BDGP6.46`.

The bundled catalogue covers every registered numbered species/release and
known dated assemblies. Some historical sources do not publish an accession;
`assembly_accession` is `None` for those entries. Unknown custom accessions can
also have unknown names or species. See the
[catalogue coverage notes](../dev/assembly-catalogue.md). Missing information is
never replaced with a guessed identifier or a mandatory network lookup.

JSON and pickle preserve the resolved metadata and cache selection. Reloading
an object does not reinterpret it against a changed catalogue. Historical
saved objects are enriched in memory without rewriting existing files.

## Discover through the new Ensembl API

Use this explicitly online entry point for a genome UUID or accession known
to the new platform:

```python
data = EnsemblRelease.from_ensembl("GCA_014441545.1")
print(data.assembly_name)    # ROS_Cfam_1.0
print(data.release)          # annotation date, currently 2020_10
print(data.genome_id)        # pinned genome UUID
print(data.website_release)  # separate website release label
```

An accession lookup selects the entry currently returned by Ensembl; a UUID
selects that specific entry. The resulting object records the concrete date,
UUID, assembly metadata, and website label for offline serialization. Species
and assembly constraints can be supplied and are checked.

Discovery uses the new metadata API's genome details and FTP links. Downloads
still use the accession-based bulk GTF/FASTA tree. It does not use the frozen
legacy REST service to find new annotations, or require GraphQL for queries.
API errors and conflicting metadata are reported rather than silently choosing
a different genome.

## Coverage, DNA, and existing caches

Existing species-based dated calls keep their defaults: human GRCh38 uses
`genes-including_alt.gtf.gz`; other assemblies use `genes.gtf.gz`. Explicit
accession selections default to primary annotation coverage. Use
`include_alt=True` only when that dataset publishes the including-alt file.
The cDNA FASTA contains all transcript biotypes. Different providers can
publish different file sets; use `Genome` for nonstandard collections.

Numbered human GRCh38 from release 82, mouse GRCm38 releases 82–102, and
zebrafish GRCz11 from release 92 use their complete `chr_patch_hapl_scaff` GTF.
Dated zebrafish has no including-alt counterpart, so its annotation coverage
is narrower. Alternative loci can add gene-name matches; prefer stable gene
IDs or a contig filter when selecting a particular locus.

Reference DNA is optional: pass `genome_fasta=True` and optionally
`genome_fasta_mask="none"`, `"soft"`, or `"hard"`. Dated datasets publish only
toplevel DNA. Read [genomic DNA](reference-dna.md) for sequence queries.

Species-based calls retain locations such as `GRCh38/ensembl2026_04/`.
Explicit-accession calls retain the original accession/provider/coverage
layout. A named assertion of the default assembly keeps the species layout;
an explicit accession uses the accession layout. Supplying
`cache_directory_path` requires a distinct directory per dataset.

Existing data and indexes are not renamed, copied, migrated, or rebuilt by
this API change. [Shared DNA caching](cache.md#how-the-shared-dna-cache-works)
continues to share new verified downloads while preserving existing private
copies. No reinstall is needed.

## Tomato numbered and dated assemblies

Tomato illustrates why numbered and dated assembly defaults are separate:
`EnsemblRelease(63, species="tomato")` selects SL4.0, while the registered dated
`EnsemblRelease("2018_10", species="tomato")` retains SL3.0 and its existing
cache path. SL4.0 publishes cDNA and peptide FASTAs without a separate ncRNA
FASTA. PyEnsembl declares that source difference explicitly; other missing
required files remain errors.

## Removed constructor

`EnsemblAnnotation` was removed in 2.30.0. Select the same dataset with:

```python
EnsemblRelease(date, assembly_accession=accession, **options)
```

The equivalent call reuses its existing cache files, including custom labels,
directories, and coverage. No reinstall or migration is needed. Imports and
serialized objects naming the removed class are no longer supported; recreate
the dataset using `EnsemblRelease`. Existing `EnsemblRelease` saved objects
continue to load.

Ensembl's [transition announcement](https://www.ensembl.info/2025/12/02/updates-to-programmatic-access-to-ensembl-and-transitioning-to-the-new-ensembl-platform/)
describes the frozen numbered services and the new platform. Both remain
accessible through this one constructor.
