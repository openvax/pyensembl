These fixtures contain unmodified Ensembl release-97 records for:

- ENSG00000285395 (XYLT1 on CHR_HG2263_PATCH), reported missing in #233.
- ENSG00000198727 (MT-CYB) and ENSG00000211459 (MT-RNR1), primary-assembly controls.

The GTF retains its metadata and all rows whose `gene_id` matches those IDs.
FASTA subsets retain complete records whose `gene:` field matches those IDs.
Coordinates and versioned identifiers are unchanged.

Sources:

- https://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.chr_patch_hapl_scaff.gtf.gz
- https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
- https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
- https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz

The full expanded GTF SHA-256 is
`8003fbf50680a944a88f7fe8dbfc67946c140ad07688b18b6def3e4ff8202675`.
The standard `Homo_sapiens.GRCh38.97.gtf.gz` contains both mitochondrial genes
but no records for ENSG00000285395.
