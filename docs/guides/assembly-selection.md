# Choose an assembly and annotation

Choose the genome assembly used by your analysis, then an annotation for that
assembly. For example, use GRCh37 annotations for positions measured on GRCh37,
not GRCh38. Annotation versions can also change gene IDs, coordinates and
sequences, so record the version you use.

| Annotation source | Choose by | Guide |
| --- | --- | --- |
| [Numbered Ensembl](https://ftp.ensembl.org/pub/) | Species, assembly and integer release | [Numbered releases](#numbered-releases) |
| [New Ensembl platform](https://www.ensembl.org/) | Species and annotation date, e.g. `2026_04` | [Dated releases](../ensembl-platform.md) |
| [Custom GTF and FASTA](custom-genomes.md) | Matched local files or URLs | [Custom genomes](custom-genomes.md) |

## Numbered releases

PyEnsembl supports numbered Ensembl releases through the final one, 116 (63 for
Ensembl Genomes); `pyensembl available` lists each species' supported range.
Ensembl's archive policy and the new platform are described in the
[platform guide](../ensembl-platform.md).
The newest supported release can change with package updates. Pin a release
number when repeating an analysis rather than relying on assembly-only selection.

```sh
pyensembl available
pyensembl install --release 93 --species human
```

`pyensembl available` lists known species, assemblies and supported release
ranges. In Python, `EnsemblRelease(93, species="human")` selects GRCh38;
`EnsemblRelease(75, species="human")` selects GRCh37. These coordinate systems
are distinct. See [Ensembl's assembly explanation](https://www.ensembl.org/info/genome/assembly/index.html).
Species accept common or Latin names, such as `"mouse"` or `"mus_musculus"`.

To choose the newest supported release for an assembly:

```sh
pyensembl install --reference-name GRCh37
```

Reference names are case-insensitive, and species is inferred from the
reference. GRCh37 selects human release 75. An explicit `--release` can select
an older compatible annotation; conflicting species, reference and release
selections are rejected. Deletion commands require an explicit release.
The Python helper `genome_for_reference_name` prefers the newest installed
release, then a downloaded release, then the newest supported release; see
[cache readiness](cache.md#list-installed-genomes).

## Annotation coverage

PyEnsembl uses Ensembl's complete `chr_patch_hapl_scaff` GTF for human GRCh38
from release 82, mouse GRCm38 releases 82–102, and zebrafish GRCz11 from release
92. These files include additional genes on assembly patches and haplotypes.
Other assemblies and earlier releases use the standard GTF filename.

Patch and haplotype contig names are preserved, for example
`CHR_HG2263_PATCH`. Gene-name searches can return additional genes on these
contigs; use stable gene IDs or a contig filter when selecting a particular locus.

If you installed one of these releases with PyEnsembl before 2.10.17, rerun
`pyensembl install` for it to add the complete annotation; existing files are
kept. Custom mirrors must provide the complete GTF filename. To use a
deliberately restricted annotation, supply its GTF as
[custom data](custom-genomes.md).

The [new platform guide](../ensembl-platform.md) explains its separate `include_alt` choice. Matching contig names alone do not prove matching assemblies.
