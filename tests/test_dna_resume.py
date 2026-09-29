"""Interrupted reference DNA downloads resume safely (#415).

A local HTTP server supports Range and If-Range with an ETag and can cut
responses short, like a dropped connection.
"""

import gzip
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import os
import random
import threading
import time

import pytest

from pyensembl.genome_fasta import GenomeFasta

# Random bases compress to ~2 MB: datacache keeps partial bytes in 1 MB
# chunks, so a cut 90% of the way through leaves bytes to resume from.
_BASES = bytes(b"ACGT"[i % 4] for i in range(256))


def fasta(seed):
    bases = random.Random(seed).randbytes(7_000_000).translate(_BASES)
    return b">1 chromosome\n" + bases + b"\n"


SEQUENCE, UPDATED_SEQUENCE = fasta(0), fasta(1)

posix_only = pytest.mark.skipif(os.name != "posix", reason="datacache resumes on POSIX only")


class RangeServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, etag):
        super().__init__(("127.0.0.1", 0), RangeHandler)
        self.etag = etag
        self.payload = gzip.compress(SEQUENCE, compresslevel=1, mtime=0)
        self.cut_short = 0  # Responses to end early, as a dropped connection.
        self.requests = []  # (Range, If-Range, status) per request.

    @property
    def url(self):
        return "http://127.0.0.1:%d/dna.fa.gz" % self.server_address[1]

    def resumed_from(self, request):
        requested, _, status = self.requests[request]
        if requested is None or status != 206:
            return 0
        return int(requested.split("=")[1].split("-")[0])


class RangeHandler(BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def do_GET(self):
        payload = self.server.payload
        requested = self.headers.get("Range")
        if_range = self.headers.get("If-Range")
        # RFC 9110: a Range with a stale If-Range gets the whole new file.
        partial = requested is not None and if_range in (None, self.server.etag)
        start = int(requested.split("=")[1].split("-")[0]) if partial else 0
        self.server.requests.append((requested, if_range, 206 if partial else 200))
        data = payload[start:]
        self.send_response(206 if partial else 200)
        if partial:
            self.send_header(
                "Content-Range", "bytes %d-%d/%d" % (start, len(payload) - 1, len(payload))
            )
        if self.server.etag:
            self.send_header("ETag", self.server.etag)
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        if self.server.cut_short:
            self.server.cut_short -= 1
            self.wfile.write(data[: len(data) * 9 // 10])
            self.close_connection = True
            return
        self.wfile.write(data)


@pytest.fixture
def server(request, monkeypatch):
    monkeypatch.setattr(time, "sleep", lambda seconds: None)  # Skip retry backoff.
    instance = RangeServer(getattr(request, "param", '"v1"'))
    thread = threading.Thread(target=instance.serve_forever, daemon=True)
    thread.start()
    yield instance
    instance.shutdown()
    instance.server_close()


def download(genome_fasta, server):
    genome_fasta._download(expected_size=len(server.payload))


def interrupt_every_attempt(genome_fasta, server):
    server.cut_short = 100
    with pytest.raises(Exception):
        download(genome_fasta, server)
    server.cut_short = 0
    return len(server.requests)


@posix_only
def test_a_dropped_connection_resumes_within_the_same_call(tmp_path, server):
    genome_fasta = GenomeFasta(server.url, tmp_path / "cache")
    server.cut_short = 1
    download(genome_fasta, server)
    assert server.resumed_from(1) > 0
    assert server.requests[1][1] == '"v1"'  # If-Range carries the ETag.
    assert genome_fasta.materialized_path.read_bytes() == SEQUENCE


@posix_only
def test_an_interrupted_download_resumes_in_a_later_call(tmp_path, server):
    genome_fasta = GenomeFasta(server.url, tmp_path / "cache")
    first_call = interrupt_every_attempt(genome_fasta, server)
    download(genome_fasta, server)
    assert server.resumed_from(first_call) > 0
    assert server.requests[first_call][1] == '"v1"'
    assert genome_fasta.materialized_path.read_bytes() == SEQUENCE
    assert not (genome_fasta.directory / "download-dna.fa.gz").exists()


@posix_only
def test_a_changed_file_restarts_instead_of_mixing_bytes(tmp_path, server):
    # Same-size uncompressed files isolate the ETag check: a new size would
    # already start a new download.
    server.payload = SEQUENCE
    genome_fasta = GenomeFasta(server.url, tmp_path / "cache")
    first_call = interrupt_every_attempt(genome_fasta, server)
    server.etag = '"v2"'
    server.payload = UPDATED_SEQUENCE
    assert len(UPDATED_SEQUENCE) == len(SEQUENCE)
    download(genome_fasta, server)
    requested, if_range, status = server.requests[first_call]
    assert requested is not None and if_range == '"v1"' and status == 200
    assert genome_fasta.materialized_path.read_bytes() == UPDATED_SEQUENCE


@posix_only
def test_a_failed_replacement_keeps_the_installed_fasta(tmp_path, server):
    genome_fasta = GenomeFasta(server.url, tmp_path / "cache")
    download(genome_fasta, server)
    interrupt_every_attempt(genome_fasta, server)
    assert genome_fasta.materialized_path.read_bytes() == SEQUENCE


@pytest.mark.parametrize("server", [None, 'W/"weak"'], indirect=True)
def test_servers_without_strong_etags_download_in_full(tmp_path, server, caplog):
    genome_fasta = GenomeFasta(server.url, tmp_path / "cache")
    download(genome_fasta, server)
    assert genome_fasta.materialized_path.read_bytes() == SEQUENCE
    if os.name == "posix":
        assert "downloading it in full" in caplog.text
