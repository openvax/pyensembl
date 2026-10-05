"""Alias lookup is explicit, ambiguity-preserving, and annotation-scoped."""

import gzip
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import threading

import datacache
import pytest
from requests.exceptions import HTTPError

from pyensembl import EnsemblRelease, GeneNameAliases, Genome


# Deliberately add a second gene to the real p53 alias to exercise ambiguity.
HGNC = (
    'symbol\talias_symbol\tprev_symbol\tensembl_gene_id\tstatus\n'
    'TP53\tp53|LFS1\tTRP53\tENSG00000141510\tApproved\n'
    'BRAF\tp53\tBRAF1\tENSG00000157764\tApproved\n'
    'MISSING\tp53\t\tENSG00000999999\tApproved\n'
    'WITHDRAWN\tp53\t\tENSG00000100000\tEntry Withdrawn\n'
    'NO_ID\tNO_GENE\t\t\tApproved\n'
)


@pytest.fixture
def genome(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        '1\ttest\tgene\t1\t10\t.\t+\t.\tgene_id "ENSG00000141510.7"; gene_name "TP53";\n'
        '1\ttest\tgene\t20\t30\t.\t+\t.\tgene_id "ENSG00000157764"; gene_name "BRAF";\n'
        '1\ttest\tgene\t40\t50\t.\t+\t.\tgene_id "AT1G01010.1"; gene_name "NAC001";\n'
    )
    with Genome("synthetic", "aliases", gtf_path_or_url=str(gtf),
                cache_directory_path=str(tmp_path / "cache")) as result:
        result.index()
        yield result


@pytest.mark.parametrize("compressed", [False, True])
def test_hgnc_alias_lookup_preserves_ambiguity_and_versions(genome, tmp_path, compressed):
    path = tmp_path / ("hgnc.tsv.gz" if compressed else "hgnc.tsv")
    path.write_bytes(gzip.compress(HGNC.encode()) if compressed else HGNC.encode())
    aliases = GeneNameAliases.from_hgnc(path)
    assert aliases.species == "homo_sapiens"
    assert aliases.source == str(path)
    assert "WITHDRAWN" not in aliases
    assert "NO_GENE" not in aliases
    assert genome.gene_ids_of_gene_name("p53", aliases=aliases) == [
        "ENSG00000141510.7", "ENSG00000157764"
    ]
    assert [gene.name for gene in genome.genes_by_name("TRP53", aliases=aliases)] == ["TP53"]
    assert genome.gene_ids_of_gene_name("TP53", aliases=aliases) == ["ENSG00000141510.7"]
    with pytest.raises(ValueError, match="Gene name not found"):
        genome.genes_by_name("P53", aliases=aliases)
    with pytest.raises(ValueError):
        genome.genes_by_name("p53")
    with pytest.raises(ValueError, match="Gene name not found"):
        genome.genes_by_name("MISSING", aliases=aliases)


def test_exact_and_alias_hits_are_combined_without_duplicates(genome):
    aliases = {"TP53": ["ENSG00000157764", "ENSG00000141510", "ENSG00000157764"]}
    assert genome.gene_ids_of_gene_name("TP53", aliases=aliases) == [
        "ENSG00000141510.7", "ENSG00000157764"
    ]
    assert genome.gene_ids_of_gene_name("TP53") == ["ENSG00000141510.7"]


def test_custom_aliases_preserve_non_ensembl_identifier_suffixes(genome):
    assert genome.gene_ids_of_gene_name("NAC", aliases={"NAC": "AT1G01010.1"}) == ["AT1G01010.1"]
    with pytest.raises(ValueError, match="Gene name not found"):
        genome.genes_by_name("NAC", aliases={"NAC": "AT1G01010"})


def test_aliases_do_not_require_gene_names_in_the_gtf(tmp_path):
    source = tmp_path / "no-names.gtf"
    source.write_text('1\ttest\tgene\t1\t10\t.\t+\t.\tgene_id "g1";\n')
    with Genome("synthetic", "no-names", gtf_path_or_url=str(source),
                cache_directory_path=str(tmp_path / "cache")) as genome:
        genome.index()
        assert genome.genes_by_name("named", aliases={"named": ["g1"]})[0].id == "g1"


def test_human_aliases_rejected_for_known_nonhuman_releases(tmp_path):
    path = tmp_path / "hgnc.tsv"
    path.write_text(HGNC)
    with pytest.raises(ValueError, match="homo_sapiens cannot be used with mus_musculus"):
        EnsemblRelease(93, species="mouse").genes_by_name(
            "p53", aliases=GeneNameAliases.from_hgnc(path)
        )


def test_mapping_copies_inputs():
    original = {"alias": ["g1"]}
    aliases = GeneNameAliases(original)
    original["alias"].append("g2")
    assert aliases["alias"] == ("g1",)
    with pytest.raises(TypeError):
        aliases["alias"] = ["g3"]


@pytest.mark.parametrize("aliases", [{"": "g1"}, {"alias": None}, {"alias": [1]}, {"alias": ""}])
def test_invalid_mapping_rejected(aliases):
    with pytest.raises(ValueError):
        GeneNameAliases(aliases)


@pytest.mark.parametrize("text", ["", "symbol\talias_symbol\n", HGNC.splitlines()[0] + "\nTP53\n",
                                  HGNC.replace("ENSG00000141510", "not-an-ensembl-id")])
def test_invalid_hgnc_data_rejected(tmp_path, text):
    path = tmp_path / "invalid.tsv"
    path.write_text(text)
    with pytest.raises(ValueError):
        GeneNameAliases.from_hgnc(path)


@pytest.fixture
def hgnc_server():
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_GET(self):
            self.server.requests.append(self.path)
            payload = self.server.payloads.get(self.path)
            if payload is None:
                self.send_error(404)
                return
            self.send_response(200)
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    server.payloads = {"/hgnc.tsv": HGNC.encode()}
    server.requests = []
    server.url = "http://127.0.0.1:%d" % server.server_address[1]
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield server
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


def test_download_hgnc_indexes_aliases_and_reuses_snapshot_offline(
    genome, hgnc_server, tmp_path, monkeypatch,
):
    source_url = hgnc_server.url + "/hgnc.tsv"
    options = dict(source_url=source_url, cache_directory_path=tmp_path / "aliases")
    aliases = GeneNameAliases.download_hgnc(**options)
    assert hgnc_server.requests == ["/hgnc.tsv"]
    assert genome.gene_ids_of_gene_name("LFS1", aliases=aliases) == ["ENSG00000141510.7"]
    assert Path(aliases.source).read_text() == HGNC
    assert datacache.inspect_file(aliases.source).source_url == source_url

    def no_network(*args, **kwargs):
        raise AssertionError("Cached aliases must not download again")

    monkeypatch.setattr(datacache, "fetch_file", no_network)
    cached = GeneNameAliases.download_hgnc(**options)
    assert cached.source == aliases.source
    assert dict(cached) == dict(aliases)


def test_download_hgnc_refresh_is_explicit(hgnc_server, tmp_path):
    options = dict(source_url=hgnc_server.url + "/hgnc.tsv", cache_directory_path=tmp_path)
    original = GeneNameAliases.download_hgnc(**options)
    hgnc_server.payloads["/hgnc.tsv"] = HGNC.replace("LFS1", "NEW_ALIAS").encode()
    cached = GeneNameAliases.download_hgnc(**options)
    assert "LFS1" in cached and "NEW_ALIAS" not in cached
    assert len(hgnc_server.requests) == 1
    refreshed = GeneNameAliases.download_hgnc(**options, overwrite=True)
    assert "NEW_ALIAS" in refreshed and "LFS1" not in refreshed
    assert "LFS1" in original
    assert len(hgnc_server.requests) == 2


def test_download_hgnc_separates_urls_with_the_same_filename(hgnc_server, tmp_path):
    hgnc_server.payloads["/archive/hgnc.tsv"] = HGNC.replace("LFS1", "OLD_ALIAS").encode()
    current = GeneNameAliases.download_hgnc(
        source_url=hgnc_server.url + "/hgnc.tsv", cache_directory_path=tmp_path,
    )
    archived = GeneNameAliases.download_hgnc(
        source_url=hgnc_server.url + "/archive/hgnc.tsv", cache_directory_path=tmp_path,
    )
    assert current.source != archived.source
    assert "LFS1" in current and "OLD_ALIAS" not in current
    assert "OLD_ALIAS" in archived and "LFS1" not in archived
    assert Path(current.source).read_text() == HGNC


@pytest.mark.parametrize("url_path", ["/hgnc.tsv.gz", "/hgnc.tsv.gz?download=1", "/snapshot"])
def test_download_hgnc_supports_gzip_snapshots(hgnc_server, tmp_path, url_path):
    hgnc_server.payloads[url_path] = gzip.compress(HGNC.encode())
    aliases = GeneNameAliases.download_hgnc(
        source_url=hgnc_server.url + url_path, cache_directory_path=tmp_path,
    )
    assert aliases["LFS1"] == ("ENSG00000141510",)
    assert Path(aliases.source).read_bytes()[:2] == b"\x1f\x8b"


def test_download_hgnc_respects_configured_cache_root(hgnc_server, tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    aliases = GeneNameAliases.download_hgnc(source_url=hgnc_server.url + "/hgnc.tsv")
    base = tmp_path / "pyensembl" / "aliases" / "homo_sapiens" / "hgnc"
    assert Path(aliases.source).is_relative_to(base)


def test_failed_hgnc_refresh_preserves_cached_snapshot(hgnc_server, tmp_path):
    options = dict(source_url=hgnc_server.url + "/hgnc.tsv", cache_directory_path=tmp_path)
    aliases = GeneNameAliases.download_hgnc(**options)
    del hgnc_server.payloads["/hgnc.tsv"]
    with pytest.raises(HTTPError):
        GeneNameAliases.download_hgnc(**options, overwrite=True)
    assert Path(aliases.source).read_text() == HGNC
    assert "LFS1" in GeneNameAliases.download_hgnc(**options)


@pytest.mark.parametrize("source_url", [None, "", "/data/hgnc.tsv"])
def test_download_hgnc_rejects_local_paths_and_missing_urls(source_url, tmp_path):
    with pytest.raises(ValueError, match="source_url must be a URL"):
        GeneNameAliases.download_hgnc(source_url=source_url, cache_directory_path=tmp_path)
    assert not list(tmp_path.iterdir())
