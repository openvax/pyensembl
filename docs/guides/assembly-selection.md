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
numbered release. These defaults have not changed. Deletion commands still
require an explicit release.

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
assemblies do not make two annotations equivalent; there is no automatic
conversion between numbered and dated releases or rolling `latest` default.
For example, human `2024_11` is a GENCODE subset without non-canonical lncRNA
transcripts; numbered release 114 contains the complete annotation.

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
