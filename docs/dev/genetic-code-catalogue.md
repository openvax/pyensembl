# Maintain genetic-code reference metadata

`pyensembl/genetic_code_catalog.json` is a pinned offline snapshot. It maps
the 43 registered Ensembl species/strain entries to verified NCBI taxon IDs
and contains 42 taxon records including their species ancestors. The wheel
and source distribution include it. Runtime construction, deserialization
and annotation never contact a taxonomy service.

## Evidence and scope

Exact species names, with accession matching when necessary, come from the
official Ensembl 116 and Ensembl Genomes 63 species TSV files. Taxonomy uses
NCBI EUtils XML: `GeneticCode/GCId`, `MitoGeneticCode/MGCId`, and plastid
`Properties/Property[PropName=pgcode]/PropValueInt`. Missing or zero values
remain unknown. Lineage records establish parent/strain relationships;
we do not treat different sibling taxa as equivalent. The platform's generic
species listing may describe one representative strain, so it is not used
to assign that strain to every reference assembly.

NCBI `gc.prt` supplies valid code IDs and names. PyEnsembl does not maintain
codon translation tables. See NCBI's [genetic-code definitions](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
and [taxonomy genetic-code documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/genetic-codes/).

The snapshot includes source URLs, SHA256 digests of uncompressed captures
and UTC retrieval time. Compressed primary inputs are checked in under
`tests/data/genetic_codes`; they allow deterministic reconstruction and
review independent of network changes. Refresh only as an intentional
metadata update. Existing serialized reference records retain their pinned
assignments after an update.

## Refresh and verify

Run from a development checkout with dependencies installed, using a new
capture directory:

```sh
python scripts/capture_genetic_code_metadata.py /tmp/ensembl-genetic-codes
python scripts/build_genetic_code_catalog.py /tmp/ensembl-genetic-codes
./lint.sh
./test.sh
```

Capture is the only networked step. It downloads bounded metadata, not
genomic files, and does not use the genome cache. The builder itself is
offline. Review species/taxon mappings, changed compartment assignments,
source digests and lineage relationships. Replace the checked-in captures
with the reviewed inputs so the reconstruction test continues to validate
the published snapshot. Additional supported species should use this source
mapping and captured NCBI evidence, not handwritten code-selection rules.
