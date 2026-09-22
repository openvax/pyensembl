"""Issue #233: include patch/haplotype annotations without reusing old indexes."""

import gzip
from pathlib import Path

import datacache
import pytest

from pyensembl import EnsemblRelease, Genome
from pyensembl.download_cache import MissingLocalFile
from pyensembl.shell import collect_selected_genomes, parser


DATA = Path(__file__).parent / "data"
PATCH_GENE = "ENSG00000285395"
PATCH_CONTIG = "CHR_HG2263_PATCH"
STANDARD_GTF = "Homo_sapiens.GRCh38.97.gtf.gz"
COMPLETE_GTF = "Homo_sapiens.GRCh38.97.chr_patch_hapl_scaff.gtf.gz"


@pytest.mark.parametrize(
    "species, release, filename",
    [
        ("human", 81, "Homo_sapiens.GRCh38.81.gtf.gz"),
        ("human", 82, "Homo_sapiens.GRCh38.82.chr_patch_hapl_scaff.gtf.gz"),
        ("human", 115, "Homo_sapiens.GRCh38.115.chr_patch_hapl_scaff.gtf.gz"),
        ("mouse", 81, "Mus_musculus.GRCm38.81.gtf.gz"),
        ("mouse", 82, "Mus_musculus.GRCm38.82.chr_patch_hapl_scaff.gtf.gz"),
        ("mouse", 102, "Mus_musculus.GRCm38.102.chr_patch_hapl_scaff.gtf.gz"),
        ("mouse", 103, "Mus_musculus.GRCm39.103.gtf.gz"),
        ("mouse", 115, "Mus_musculus.GRCm39.115.gtf.gz"),
        ("danio_rerio", 91, "Danio_rerio.GRCz10.91.gtf.gz"),
        ("danio_rerio", 92, "Danio_rerio.GRCz11.92.chr_patch_hapl_scaff.gtf.gz"),
        ("danio_rerio", 115, "Danio_rerio.GRCz11.115.chr_patch_hapl_scaff.gtf.gz"),
        ("rat", 97, "Rattus_norvegicus.Rnor_6.0.97.gtf.gz"),
        ("arabidopsis_thaliana", 58, "Arabidopsis_thaliana.TAIR10.58.gtf.gz"),
    ],
)
def test_gtf_selection_matches_archived_ensembl_files(species, release, filename):
    # These filenames were verified in the official archive, including the
    # mouse switch back to standard filenames at the GRCm39 transition.
    genome = EnsemblRelease(release, species=species)
    assert genome.gtf_url.endswith("/" + filename)
    assert "/release-%d/" % release in genome.gtf_url


def write_compressed(path, text):
    with gzip.open(path, "wt") as output:
        output.write(text)


def test_mirror_missing_complete_gtf_does_not_fall_back(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "cache"))
    write_compressed(tmp_path / STANDARD_GTF, "# legacy annotation\n")
    genome, = collect_selected_genomes(parser.parse_args([
        "install", "--release", "97", "--custom-mirror", str(tmp_path),
    ]))
    with pytest.raises(MissingLocalFile) as error:
        genome.download()
    assert Path(error.value.path).name == COMPLETE_GTF


@pytest.mark.parametrize("custom_mirror", [False, True])
def test_reported_patch_gene_survives_selection_indexing_and_cached_upgrade(
    tmp_path, monkeypatch, custom_mirror,
):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "cache"))

    def no_network(*args, **kwargs):
        pytest.fail("The source-grounded regression fixture must work offline")

    monkeypatch.setattr(datacache.download, "_download_and_decompress_if_necessary", no_network)
    release = EnsemblRelease(97)
    cache = Path(release.download_cache.cache_directory_path)
    cache.mkdir(parents=True)
    sources = tmp_path / "mirror" if custom_mirror else cache
    sources.mkdir(exist_ok=True)
    complete_gtf = (DATA / "human.ensembl.97.patch.gtf").read_text()
    # The standard upstream GTF contains the mitochondrial controls but none
    # of the reported patch gene's rows. Keep both files, as on an upgrade.
    standard_gtf = "".join(
        line for line in complete_gtf.splitlines(keepends=True)
        if PATCH_GENE not in line
    )
    write_compressed(sources / STANDARD_GTF, standard_gtf)
    write_compressed(sources / COMPLETE_GTF, complete_gtf)
    for suffix, fixture in [
        ("cdna.all.fa.gz", "cdna.fa"),
        ("ncrna.fa.gz", "ncrna.fa"),
        ("pep.all.fa.gz", "pep.fa"),
    ]:
        write_compressed(
            sources / ("Homo_sapiens.GRCh38." + suffix),
            (DATA / ("human.ensembl.97.patch." + fixture)).read_text(),
        )

    with Genome(
        "GRCh38", "ensembl", 97,
        gtf_path_or_url=str(sources / STANDARD_GTF),
        cache_directory_path=str(cache),
    ) as previous:
        previous.index()
        assert PATCH_GENE not in previous.gene_ids()
        old_index = previous.db.local_db_path

    arguments = ["install", "--release", "97"]
    if custom_mirror:
        arguments += ["--custom-mirror", str(sources)]
    genome, = collect_selected_genomes(parser.parse_args(arguments))
    with genome:
        genome.download()
        genome.index()
        assert {gene.gene_id for gene in genome.genes()} == {
            "ENSG00000198727", "ENSG00000211459", PATCH_GENE,
        }
        gene = genome.gene_by_id(PATCH_GENE)
        assert gene.gene_name == "XYLT1"
        assert (gene.contig, gene.start, gene.end, gene.strand) == (
            PATCH_CONTIG, 17101769, 17470881, "-",
        )
        assert len(gene.transcripts) == 4
        assert genome.gene_ids_at_locus(PATCH_CONTIG, 17470796) == [PATCH_GENE]
        assert set(genome.contigs()) == {"MT", PATCH_CONTIG}
        transcript = genome.transcript_by_id("ENST00000644858")
        assert len(transcript.exons) == 12
        assert len(transcript.sequence) == sum(len(exon) for exon in transcript.exons)
        assert transcript.protein_id == "ENSP00000495604"
        assert transcript.protein_sequence
        assert genome.db.local_db_path != old_index
        assert Path(old_index).exists()
