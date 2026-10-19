"""A CLI mirror install can be reopened without acquisition in another process."""
from functools import partial
import gzip
import hashlib
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import os
from pathlib import Path
import pickle
import subprocess
import sys
from threading import Thread

import pytest

from pyensembl import EnsemblRelease, Genome
from pyensembl.shell import all_combinations_of_ensembl_genomes, parser


def snapshot(directory):
    return {str(p.relative_to(directory)): (hashlib.sha256(p.read_bytes()).hexdigest(),
                                          p.stat().st_mtime_ns, p.stat().st_ctime_ns,
                                          p.stat().st_ino, p.stat().st_mode)
            for p in directory.rglob("*") if p.is_file()}


def test_cli_install_reopens_from_python_in_fresh_process(tmp_path):
    source = tmp_path / "mirror"
    source.mkdir()
    fixtures = Path(__file__).parent / "data"
    files = {
        "Mus_musculus.GRCm38.81.gtf.gz": "mouse.ensembl.81.partial.ENSMUSG00000017167.gtf",
        "Mus_musculus.GRCm38.cdna.all.fa.gz": "mouse.ensembl.81.partial.ENSMUSG00000017167.fa",
        "Mus_musculus.GRCm38.ncrna.fa.gz": "mouse.ensembl.81.partial.ncrna.ENSMUSG00000017167.fa",
        "Mus_musculus.GRCm38.pep.all.fa.gz": "mouse.ensembl.81.partial.ENSMUSG00000017167.pep",
    }
    for name, fixture in files.items():
        (source / name).write_bytes(gzip.compress((fixtures / fixture).read_bytes()))
    server = ThreadingHTTPServer(("127.0.0.1", 0), partial(SimpleHTTPRequestHandler, directory=str(source)))
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    mirror = "http://127.0.0.1:%d/" % server.server_port
    cache = tmp_path / "cache"
    env = {**os.environ, "PYENSEMBL_CACHE_DIR": str(cache)}
    try:
        installed = subprocess.run(
            [sys.executable, "-c", "from pyensembl.shell import run; run()",
             "install", "--release", "81", "--species", "mouse", "--custom-mirror", mirror],
            env=env, capture_output=True, text=True,
        )
        assert installed.returncode == 0, installed.stderr
    finally:
        server.shutdown()
        server.server_close()
        thread.join()
    before = snapshot(cache)
    assert any(name.endswith(".db") for name in before)
    reopened = subprocess.run([sys.executable, "-c", """
import sys
from pyensembl import EnsemblRelease, Genome
from pyensembl.download_cache import DownloadCache
from pyensembl.database import Database
import pyensembl.sequence_data as sequence_data

def forbid(*args, **kwargs):
    raise AssertionError('Reopening installed mirror data must not download or index')

DownloadCache._fetch = forbid
Database.create = forbid
sequence_data.parse_fasta_dictionary = forbid
Genome.index = forbid
canonical = EnsemblRelease(81, species='mouse')
genome = canonical.with_mirror(sys.argv[1])
assert genome.installed()
assert genome.gene_by_id('ENSMUSG00000017167').name == 'Cntnap1'
assert genome.transcript_by_id('ENSMUST00000103109').sequence
assert genome.transcript_by_id('ENSMUST00000103109').protein_sequence.startswith('MMSLRL')
assert not canonical.installed()
assert not canonical.with_mirror(sys.argv[1] + 'another/').installed()
restored = Genome.from_json(genome.to_json())
assert restored.installed()
assert restored.transcript_by_id('ENSMUST00000103109').protein_sequence.startswith('MMSLRL')
""", mirror], env=env, capture_output=True, text=True)
    assert reopened.returncode == 0, reopened.stderr
    assert snapshot(cache) == before


@pytest.mark.parametrize("release", [75, "2026_04"])
def test_cli_and_python_use_same_mirror_sources_and_roundtrip(release):
    mirror = "https://mirror.example/release/"
    api = EnsemblRelease(release, genome_fasta=True).with_mirror(mirror)
    cli, = all_combinations_of_ensembl_genomes(parser.parse_args([
        "install", "--release", str(release), "--custom-mirror", mirror,
        "--with-genome-fasta",
    ]))
    assert cli == api
    assert isinstance(api, Genome)
    assert not isinstance(api, EnsemblRelease)
    assert api._genome_fasta_path_or_url.startswith(mirror)
    assert Genome.from_json(api.to_json()) == api
    assert pickle.loads(pickle.dumps(api)) == api


def test_mirror_preserves_explicit_dna_source(tmp_path):
    local = tmp_path / "dna.fa"
    local.write_text(">1\nACGT\n")
    mirror = "https://mirror.example/"
    genome = EnsemblRelease(75, genome_fasta=local).with_mirror(mirror)
    assert genome._genome_fasta_path_or_url == str(local)
    remote = "https://dna.example/reference.fa.gz"
    override = EnsemblRelease(75, genome_fasta=True).with_mirror(
        mirror, genome_fasta_path_or_url=remote)
    assert override._genome_fasta_path_or_url == remote


def test_empty_mirror_rejected():
    with pytest.raises(ValueError, match="nonempty"):
        EnsemblRelease(75).with_mirror("")
