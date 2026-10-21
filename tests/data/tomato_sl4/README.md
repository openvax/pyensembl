# Tomato source and compatibility fixtures

`genes.gtf`, `cdna.fa` and `pep.fa` are verbatim subsets of Ensembl Genomes 63
tomato SL4.0, retrieved on 2026-10-20. The two transcripts are
`mRNA-Solyc04g011850.1.1` and `mRNA-Solyc02g090330.3.1`. `provenance.json`
records the upstream GTF URL and SHA-256 digests of the complete compressed
GTF/cDNA/peptide inputs. FASTAs are under the corresponding
`/pub/release-63/plants/fasta/solanum_lycopersicum_gca000188115v5cm/{cdna,pep}/`
directories. The independent GTF header capture is in `../assembly_metadata/`.

`pyensembl-2.30.0-cache.tar.gz` and its JSON describe real tiny GTF/FASTA/SQLite
installations created with unmodified 2.30.0 release code, before this change.
Their synthetic gene/transcript/protein markers identify old numbered SL3.0
release 62 and dated SL3.0 `2018_10` independently. Only HTTP transport was
replaced; acquisition, parsing and indexing used the released implementation.
Do not regenerate these compatibility fixtures with newer PyEnsembl code.
