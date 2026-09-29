"""An interrupted reference DNA download resumes in a later attempt."""

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
SEQUENCE = b">1 chromosome\n" + random.Random(0).randbytes(7_000_000).translate(_BASES) + b"\n"
PAYLOAD = gzip.compress(SEQUENCE, compresslevel=1, mtime=0)


class RangeServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, etag):
        super().__init__(("127.0.0.1", 0), RangeHandler)
        self.etag = etag
        self.cut_short = 0  # Responses to end early, as a dropped connection.
        self.ranges = []


class RangeHandler(BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def do_GET(self):
        requested = self.headers.get("Range")
        self.server.ranges.append(requested)
        start = int(requested.split("=")[1].split("-")[0]) if requested else 0
        data = PAYLOAD[start:]
        self.send_response(206 if requested else 200)
        if requested:
            self.send_header(
                "Content-Range", "bytes %d-%d/%d" % (start, len(PAYLOAD) - 1, len(PAYLOAD))
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


def url(server):
    return "http://127.0.0.1:%d/dna.fa.gz" % server.server_address[1]


@pytest.mark.skipif(os.name != "posix", reason="datacache resumes on POSIX only")
def test_interrupted_download_resumes_in_a_later_attempt(tmp_path, server):
    genome_fasta = GenomeFasta(url(server), tmp_path / "cache")
    server.cut_short = 100  # Every retry of the first attempt is cut off.
    with pytest.raises(Exception):
        genome_fasta._download(expected_size=len(PAYLOAD))
    server.cut_short = 0
    first_attempt = len(server.ranges)
    genome_fasta._download(expected_size=len(PAYLOAD))
    resumed_from = server.ranges[first_attempt]
    assert resumed_from is not None and int(resumed_from.split("=")[1].split("-")[0]) > 0
    assert genome_fasta.materialized_path.read_bytes() == SEQUENCE
    assert not (genome_fasta.directory / "download-dna.fa.gz").exists()


@pytest.mark.parametrize("server", [None, 'W/"weak"'], indirect=True)
def test_servers_without_strong_etags_download_in_full(tmp_path, server, caplog):
    genome_fasta = GenomeFasta(url(server), tmp_path / "cache")
    genome_fasta._download(expected_size=len(PAYLOAD))
    assert genome_fasta.materialized_path.read_bytes() == SEQUENCE
    if os.name == "posix":
        assert "downloading it in full" in caplog.text
