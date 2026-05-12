"""
Regression test for issue #219: when a Genome is constructed with a local GTF
and no explicit cache_directory_path, the sqlite database should land in the
download cache rather than next to the GTF file (which may be read-only).
"""

import os

from pyensembl import Genome

from .common import TemporaryDirectory


GTF_BODY = (
    '1\ttest\ttranscript\t100\t500\t.\t+\t.\t'
    'gene_id "G1"; transcript_id "T1"; gene_name "FN1";\n'
    '1\ttest\texon\t100\t200\t.\t+\t.\t'
    'gene_id "G1"; transcript_id "T1"; exon_number "1"; gene_name "FN1";\n'
    '1\ttest\texon\t300\t500\t.\t+\t.\t'
    'gene_id "G1"; transcript_id "T1"; exon_number "2"; gene_name "FN1";\n'
)


def test_genome_db_cache_falls_back_to_download_cache(monkeypatch):
    """Without an explicit cache_directory_path the Database should be
    written under the download cache, not alongside the source GTF."""
    with TemporaryDirectory() as cache_root, TemporaryDirectory() as gtf_dir:
        monkeypatch.setenv("PYENSEMBL_CACHE_DIR", cache_root)
        gtf_path = os.path.join(gtf_dir, "fix_219.gtf")
        with open(gtf_path, "w") as f:
            f.write(GTF_BODY)

        genome = Genome(
            reference_name="GRCh38",
            annotation_name="fix_219_test",
            gtf_path_or_url=gtf_path,
        )

        download_cache_dir = genome.download_cache.cache_directory_path
        assert download_cache_dir.startswith(cache_root)
        assert genome.db.cache_directory_path == download_cache_dir
        assert genome.db.cache_directory_path != gtf_dir


def test_genome_db_respects_explicit_cache_directory():
    """An explicit cache_directory_path is still honored end-to-end."""
    with TemporaryDirectory() as cache_dir, TemporaryDirectory() as gtf_dir:
        gtf_path = os.path.join(gtf_dir, "fix_219_explicit.gtf")
        with open(gtf_path, "w") as f:
            f.write(GTF_BODY)

        genome = Genome(
            reference_name="GRCh38",
            annotation_name="fix_219_test_explicit",
            gtf_path_or_url=gtf_path,
            cache_directory_path=cache_dir,
        )
        assert genome.db.cache_directory_path == cache_dir
        assert genome.db.cache_directory_path != gtf_dir
