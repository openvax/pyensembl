"""Install the captured upstream SL4.0 files, and retain old SL3.0 caches."""

import gzip
from io import BytesIO
import json
from pathlib import Path
import pickle
import shlex
import tarfile

import pytest
import requests

from pyensembl import EnsemblRelease, shell
from pyensembl.species import Species, tomato
from .test_legacy_cache_compatibility import snapshot
from .test_release_cache_combinations import forbid_parsing


DATA = Path(__file__).parent / "data" / "tomato_sl4"
ACCESSION = "GCA_000188115.5"


@pytest.fixture
def upstream(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "cache"))
    calls = []

    def request(session, method, url, **kwargs):
        assert method.upper() == "GET"
        assert "/ncrna/" not in url
        name = "genes.gtf" if "/gtf/" in url else "cdna.fa" if "/cdna/" in url else "pep.fa"
        payload = gzip.compress((DATA / name).read_bytes(), mtime=0)
        calls.append(url)
        response = requests.Response()
        response.status_code = 200
        response.url = url
        response.headers["Content-Length"] = str(len(payload))
        response.raw = BytesIO(payload)
        return response

    monkeypatch.setattr(requests.sessions.Session, "request", request)
    return calls


@pytest.mark.parametrize("options", [
    {"species": "tomato"}, {"species": "solanum_lycopersicum"},
    {"species": "solanum_lycopersicum_gca000188115v5cm"},
    {"assembly": "SL4.0"}, {"assembly_accession": ACCESSION},
])
def test_install_and_reopen_real_tomato_slice(upstream, tmp_path, monkeypatch, options):
    with EnsemblRelease(63, **options) as genome:
        assert (genome.assembly, genome.assembly_name, genome.assembly_accession) == (
            "SL4.0", "SL4.0", ACCESSION,
        )
        genome.download()
        genome.index()
        assert genome.installed()
        assert len(genome.gene_ids()) == 2
        assert len(genome.transcript_by_id("mRNA-Solyc04g011850.1.1").sequence) == 309
        assert len(genome.transcript_by_id("mRNA-Solyc02g090330.3.1").protein_sequence) == 228
        encoded, pickled, command = genome.to_json(), pickle.dumps(genome), genome.install_string()
        paths = genome.gtf_path, genome.db.local_db_path, genome.transcript_fasta_paths
    assert len(upstream) == 3
    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **k: pytest.fail("network"))
    forbid_parsing(monkeypatch)
    root = tmp_path / "cache"
    before = snapshot(root)
    for restored in [
        EnsemblRelease(63, species="tomato"), EnsemblRelease.from_json(encoded),
        pickle.loads(pickled),
        shell.collect_selected_genomes(shell.parser.parse_args(shlex.split(command)[1:]))[0],
    ]:
        with restored:
            restored.download()
            restored.index()
            assert (restored.gtf_path, restored.db.local_db_path, restored.transcript_fasta_paths) == paths
            assert restored.transcript_by_id("mRNA-Solyc04g011850.1.1").sequence.startswith("ATGGATATG")
    assert snapshot(root) == before


def test_absent_required_cdna_is_not_silently_ignored(upstream, monkeypatch):
    request = requests.sessions.Session.request

    def missing(session, method, url, **kwargs):
        if "/cdna/" in url:
            raise FileNotFoundError("required cDNA is missing")
        return request(session, method, url, **kwargs)

    monkeypatch.setattr(requests.sessions.Session, "request", missing)
    with EnsemblRelease(63, species="tomato") as genome:
        with pytest.raises(Exception, match="required cDNA is missing"):
            genome.download()
        assert not genome.installed()


def test_numbered_and_dated_tomato_remain_distinct():
    old = EnsemblRelease(62, species="tomato")
    dated = EnsemblRelease("2018_10", species="tomato")
    current = EnsemblRelease(63, species="tomato")
    assert old.reference_name == dated.reference_name == "SL3.0"
    assert dated.assembly_accession == "GCA_000188115.3"
    assert current.reference_name == "SL4.0"
    assert any("/ncrna/" in url for url in old.transcript_fasta_urls)
    assert "/solanum_lycopersicum/" in old.gtf_url
    assert len(current.transcript_fasta_urls) == 1
    mirror = current.with_mirror("https://mirror.example/ensembl")
    assert len(mirror._transcript_fasta_paths_or_urls) == 1
    assert "/ncrna/" not in mirror._transcript_fasta_paths_or_urls[0]
    assert tomato.which_reference("2018_10") == "SL3.0"
    table = shell.format_available_species(use_color=False, dated_releases={tomato.latin_name: ["2018_10"]})
    assert "SL3.0" in next(line for line in table.splitlines() if "2018_10" in line)


def test_source_declaration_distinguishes_cached_custom_species(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    options = dict(latin_name="test_source", reference_assemblies={"Test": (63, 63)})
    ordinary = Species(**options)
    without = Species(**options, assemblies_without_ncrna=("Test",))
    assert ordinary != without
    assert len({ordinary: 1, without: 2}) == 2
    first = EnsemblRelease.cached(63, species=ordinary)
    second = EnsemblRelease.cached(63, species=without)
    assert first is not second
    assert len(first.transcript_fasta_urls) == 2
    assert len(second.transcript_fasta_urls) == 1


def test_published_old_tomato_caches_reopen_read_only(tmp_path, monkeypatch):
    with tarfile.open(DATA / "pyensembl-2.30.0-cache.tar.gz") as archive:
        archive.extractall(tmp_path, filter="data")
    root = tmp_path / "cache"
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(root))
    monkeypatch.setattr(requests.sessions.Session, "request", lambda *a, **kw: pytest.fail("network"))
    forbid_parsing(monkeypatch)
    modes = {path: path.stat().st_mode for path in root.rglob("*")}
    for path in modes:
        path.chmod(0o555 if path.is_dir() else 0o444)
    before = snapshot(root)
    try:
        for saved in json.loads((DATA / "pyensembl-2.30.0.json").read_text()):
            for genome in [EnsemblRelease(saved["release"], species="tomato"),
                           EnsemblRelease.from_json(saved["json"])]:
                with genome:
                    genome.download()
                    genome.index()
                    assert genome.transcript_by_id(saved["expected"]["transcript"]).sequence == saved["expected"]["cdna"]
                    assert genome.reference_name == "SL3.0"
        assert any(g.species == tomato and g.release == "2018_10" for g in shell._cached_dated_releases())
        assert snapshot(root) == before
    finally:
        for path, mode in modes.items():
            path.chmod(mode)
