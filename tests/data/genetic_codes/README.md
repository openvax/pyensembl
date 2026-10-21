# Genetic-code metadata captures

The compressed TSV, XML and gc.prt files are verbatim upstream responses,
gzip-compressed with mtime=0. `capture.json` records the UTC retrieval time
and complete NCBI request. The generated package catalogue records source
URLs and SHA256 digests of the uncompressed responses.

The Ensembl 116 and Ensembl Genomes 63 species files map exact registered
reference names/accessions to source taxon IDs. NCBI Taxonomy records include
those taxa and their species ancestors. In particular, yeast S288C and
Candida SC5314 retain their strain taxon IDs; domestic dog retains taxon
9615 and its species ancestor 9612. The shared Ensembl platform species
listing's representative strain is not used as a reference taxon.

`scripts/build_genetic_code_catalog.py` builds from these captures offline.
`scripts/capture_genetic_code_metadata.py` is the separate, explicitly
networked maintenance command. Capture again only when intentionally updating
and reviewing the taxonomy snapshot.

`pyensembl-2.30.1-cache.tar.gz` and its JSON settings were created using the
unmodified published 2.30.1 implementation. Its tiny synthetic GTF declares
table 5 on a CDS, while its original SQLite index omits that attribute. This
fixture proves the new code can retrieve a formerly unindexed override
without rewriting an old cache. It contains no biological sequence data.
