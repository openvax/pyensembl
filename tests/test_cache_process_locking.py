"""Real process contention and recovery for shared DNA cache locks."""

from contextlib import contextmanager
import gzip
import multiprocessing
import os
from pathlib import Path
from unittest.mock import patch

from filelock import Timeout
import pytest

from pyensembl import EnsemblRelease
from pyensembl import genome_fasta, genome_fasta_cache as cache


PAYLOAD = gzip.compress(b">1\nACGTACGT\n")
IDENTITY = {
    "provider": "ftp.ensembl.org", "species": "homo_sapiens",
    "assembly": "GCA_000001405.18", "filename": "Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
    "checksum": [1234, 1], "compressed_size": len(PAYLOAD),
}


@contextmanager
def controlled_download(root, entered=None, finish=None):
    def fetch(source, destination, **kwargs):
        with (Path(root) / "downloads.txt").open("a") as log:
            log.write(source + "\n")
        if entered is not None:
            Path(destination).write_bytes(PAYLOAD[:5])
            (Path(destination).parent / ".interrupted-staging").write_bytes(b"partial")
            entered.set()
            if not finish.wait(30):
                raise RuntimeError("Timed out waiting for the test coordinator")
        Path(destination).write_bytes(PAYLOAD)
        return str(destination)

    with patch.dict(os.environ, {"PYENSEMBL_CACHE_DIR": str(root)}), \
            patch.object(cache, "_remote_identity", lambda source: dict(IDENTITY)), \
            patch.object(genome_fasta, "fetch_file", fetch):
        yield


def install_worker(root, release, start, entered, finish, results):
    try:
        with controlled_download(root, entered, finish):
            results.put("ready")
            if not start.wait(30):
                raise RuntimeError("Test did not start installation")
            with EnsemblRelease(release, genome_fasta=True) as genome:
                genome.download_genome_fasta()
                genome.index_genome_fasta()
                results.put((genome.genome_fasta_path, genome.sequence("1", 1, 8)))
    except Exception as error:
        results.put(repr(error))
        raise


def stop_workers(workers):
    for worker in workers:
        if worker.is_alive():
            worker.terminate()
        worker.join(10)
        assert not worker.is_alive()


def test_simultaneous_download_and_index_callers_share_one_object(tmp_path):
    context = multiprocessing.get_context("spawn")
    start, entered, finish = (context.Event() for _ in range(3))
    results = context.Queue()
    workers = [context.Process(target=install_worker, args=(
        str(tmp_path), release, start, entered, finish, results
    )) for release in (81, 82)]
    try:
        for worker in workers:
            worker.start()
        assert [results.get(timeout=30) for _ in workers] == ["ready", "ready"]
        start.set()
        assert entered.wait(30)
        root = tmp_path / "pyensembl" / "dna_cache"
        directory = cache._object_directory(root, IDENTITY)
        with pytest.raises(Timeout):
            with cache._object_lock(root, directory, timeout=0.1):
                pytest.fail("A second process acquired an active writer's lock")
        finish.set()
        values = [results.get(timeout=30) for _ in workers]
        assert values[0] == values[1]
        assert values[0][1] == "ACGTACGT"
        for worker in workers:
            worker.join(30)
            assert worker.exitcode == 0
        assert len((tmp_path / "downloads.txt").read_text().splitlines()) == 1
    finally:
        finish.set()
        stop_workers(workers)
        results.close()


def test_terminated_writer_releases_lock_and_cache_can_be_reused(tmp_path):
    context = multiprocessing.get_context("spawn")
    start, entered, finish = (context.Event() for _ in range(3))
    results = context.Queue()
    worker = context.Process(target=install_worker, args=(
        str(tmp_path), 81, start, entered, finish, results
    ))
    try:
        worker.start()
        assert results.get(timeout=30) == "ready"
        start.set()
        assert entered.wait(30)
        worker.terminate()
        worker.join(10)
        assert not worker.is_alive()
        with controlled_download(tmp_path):
            with EnsemblRelease(81, genome_fasta=True) as genome:
                assert genome.genome_fasta_path is None
                genome.download_genome_fasta()
                genome.index_genome_fasta()
                assert genome.sequence("1", 1, 8) == "ACGTACGT"
                assert not (genome._genome_fasta.directory / ".interrupted-staging").exists()
                # Ready data can be reopened without downloading or acquiring
                # a writer lock; an interrupted installer leaves no live hold.
                with patch.object(genome_fasta, "fetch_file", side_effect=AssertionError):
                    with EnsemblRelease(81, genome_fasta=True) as reopened:
                        reopened.download_genome_fasta()
                        assert reopened.sequence("1", 1, 8) == "ACGTACGT"
    finally:
        stop_workers([worker])
        results.close()
