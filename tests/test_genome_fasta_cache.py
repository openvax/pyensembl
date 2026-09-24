"""Shared DNA cache tests; real files/indexes with controlled HTTP metadata."""

import gzip
import io
import json
import os
from pathlib import Path
import shutil
import stat
import threading

import pytest

from pyensembl import EnsemblRelease, Genome, MissingGenomeFastaError, prune_genome_fastas
from pyensembl import genome_fasta_cache as cache
from pyensembl.genome_fasta import GenomeFasta
from .test_genome_fasta import DNA, forbid_downloads, run_cli, serve_downloads

# README preambles as published: Ensembl labels the accession, while Ensembl
# Genomes puts it alone on the next line.
ENSEMBL_README = (
    "#### README ####\n\nThe genome assembly represented here corresponds to "
    "GenBank Assembly ID \n%s\n\n#######################\nFasta DNA dumps\n"
)
ENSEMBL_GENOMES_README = (
    "#### README ####\n\nThe genome assembly represented here corresponds to  \n"
    "%s\n\n#######################\nFasta DNA dumps\n"
)


@pytest.fixture
def shared_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    calls = serve_downloads(monkeypatch, tmp_path, lambda url: gzip.compress(DNA))

    def metadata(request, **kwargs):
        url = request if isinstance(request, str) else request.full_url
        release = int(url.split("release-")[1].split("/")[0])
        if url.endswith("README"):
            if "ensemblgenomes" in url:
                return io.BytesIO((ENSEMBL_GENOMES_README % "GCA_000001735.1").encode())
            assembly = "GCA_000001405.18" if release <= 82 else "GCA_000001405.20"
            return io.BytesIO((ENSEMBL_README % assembly).encode())
        if url.endswith("CHECKSUMS"):
            lines = [
                "17918 1010762 %s.%s.%s.fa.gz" % (prefix, mask, flavor)
                for prefix in ("Homo_sapiens.GRCh38", "Arabidopsis_thaliana.TAIR10")
                for mask in ("dna", "dna_sm", "dna_rm")
                for flavor in ("toplevel", "primary_assembly")
            ]
            return io.BytesIO("\n".join(lines).encode())
        assert request.get_method() == "HEAD"
        response = io.BytesIO()
        response.headers = {"Content-Length": str(len(gzip.compress(DNA)))}
        return response

    monkeypatch.setattr(cache, "urlopen", metadata)
    return calls


def release(number=81, **kwargs):
    return EnsemblRelease(number, download_genome_fasta=True, **kwargs)


def test_same_patch_reuses_file_and_index_then_works_offline(shared_cache, monkeypatch):
    first, second = release(81), release(82)
    first.download_genome_fasta()
    first.index_genome_fasta()
    first_index = first._genome_fasta.index_path.stat()
    second.download_genome_fasta()
    second.index_genome_fasta()
    assert len(shared_cache) == 1
    assert first.genome_fasta_path == second.genome_fasta_path
    assert first._genome_fasta.index_path == second._genome_fasta.index_path
    assert second._genome_fasta.index_path.stat().st_mtime_ns == first_index.st_mtime_ns
    assert first.sequence("MT", 1, 4) == second.sequence("MT", 1, 4) == "GCTA"
    assert first._genome_fasta.identity["assembly"] == "GCA_000001405.18"
    relative = Path(first.genome_fasta_path).relative_to(cache.dna_cache_root())
    assert relative.parts[:-2] == (
        "homo_sapiens",
        "ftp.ensembl.org",
        "GRCh38-GCA_000001405.18",
        "toplevel",
        "unmasked",
        "fasta",
    )
    assert len(relative.parts[-2]) == 16
    assert relative.name == "sequence.fa"
    first.close()
    second.close()

    def offline(*args, **kwargs):
        pytest.fail("Installed shared DNA should work offline")

    monkeypatch.setattr(cache, "urlopen", offline)
    forbid_downloads(monkeypatch, "Installed shared DNA should work offline")
    restored = release(82)
    restored.download_genome_fasta()
    assert restored.sequence("1", 3, 10) == "GTNNTTAA"
    restored.close()


def test_patch_mask_and_flavor_separate_objects(shared_cache):
    genomes = [
        release(81),
        release(83),
        release(81, genome_fasta_mask="soft"),
        release(81, genome_fasta_type="primary_assembly"),
    ]
    for genome in genomes:
        genome.download_genome_fasta()
        assert genome.sequence("MT", 1, 4) == "GCTA"
        genome.close()
    assert len(shared_cache) == 4
    assert len({genome.genome_fasta_path for genome in genomes}) == 4
    assert "GRCh38-GCA_000001405.20" in Path(genomes[1].genome_fasta_path).parts
    assert "softmasked" in Path(genomes[2].genome_fasta_path).parts
    assert "primary_assembly" in Path(genomes[3].genome_fasta_path).parts
    # All installed flavors remain referenced, not just the most recent one.
    assert prune_genome_fastas() == []
    assert Path(release(81).genome_fasta_path).exists()


def test_changed_upstream_file_metadata_separates_same_patch(shared_cache, monkeypatch):
    first = release(81)
    first.download_genome_fasta()
    original = cache.urlopen

    def changed(request, **kwargs):
        response = original(request, **kwargs)
        if isinstance(request, str) and request.endswith("CHECKSUMS"):
            response = io.BytesIO(response.read().replace(b"17918", b"17919"))
        return response

    monkeypatch.setattr(cache, "urlopen", changed)
    second = release(82)
    second.download_genome_fasta()
    assert first.genome_fasta_path != second.genome_fasta_path
    assert (
        Path(first.genome_fasta_path).parent.parent
        == Path(second.genome_fasta_path).parent.parent
    )
    assert len(shared_cache) == 2


def test_incomplete_metadata_uses_release_specific_cache(shared_cache, monkeypatch):
    monkeypatch.setattr(
        cache, "urlopen", lambda *a, **k: io.BytesIO(b"no assembly metadata")
    )
    first, second = release(81), release(82)
    first.download_genome_fasta()
    second.download_genome_fasta()
    assert first.genome_fasta_path != second.genome_fasta_path
    assert "GRCh38-unverified" in Path(first.genome_fasta_path).parts
    assert len(shared_cache) == 2


def test_local_and_mirror_sources_do_not_enter_shared_cache(shared_cache, tmp_path):
    path = tmp_path / "custom.fa"
    path.write_bytes(DNA)
    local = release(81, genome_fasta_path=path)
    mirror = release(81, server="https://example.test")
    for genome in (local, mirror):
        genome.download_genome_fasta()
        genome.index_genome_fasta()
        assert not isinstance(genome._genome_fasta, cache.SharedGenomeFasta)
        assert genome.sequence("MT", 1, 4) == "GCTA"
        genome.close()
    assert not cache.dna_cache_root().exists()
    assert prune_genome_fastas() == []
    assert path.read_bytes() == DNA


def test_prune_preserves_live_references_and_supports_preview(
    shared_cache, monkeypatch, capsys
):
    first, second = release(81), release(82)
    for genome in (first, second):
        genome.download_genome_fasta()
        genome.index_genome_fasta()
        genome.close()
    object_path = first._genome_fasta.directory
    run_cli(monkeypatch, "delete-index-files", "--release", "81")
    assert first._genome_fasta.index_path.exists()  # Shared index still serves 82.
    run_cli(monkeypatch, "delete-all-files", "--release", "81")
    assert prune_genome_fastas() == []
    assert second.sequence("MT", 1, 4) == "GCTA"
    second.close()
    run_cli(monkeypatch, "delete-all-files", "--release", "82")
    preview = prune_genome_fastas(dry_run=True)
    assert len(preview) == 1 and preview[0][0] == str(object_path) and preview[0][1] > 0
    assert object_path.exists()
    run_cli(monkeypatch, "prune", "--orphan-genome-fastas", "--dry-run")
    assert "Would delete" in capsys.readouterr().out
    assert object_path.exists()
    run_cli(monkeypatch, "prune", "--orphan-genome-fastas")
    assert "Deleted" in capsys.readouterr().out
    assert not object_path.exists()
    assert prune_genome_fastas() == []


def test_malformed_references_abort_pruning_before_any_deletion(shared_cache):
    genome = release()
    genome.download_genome_fasta()
    manifest = genome._genome_fasta.reference_path
    manifest.write_text("not json")
    with pytest.raises(ValueError, match="Cannot safely prune"):
        prune_genome_fastas()
    assert Path(genome.genome_fasta_path).exists()
    manifest.write_text(
        json.dumps({"source": "url", "shared_key": "../outside", "identity": {}})
    )
    with pytest.raises(ValueError, match="invalid shared reference"):
        prune_genome_fastas()
    with pytest.raises(ValueError, match="Invalid shared"):
        # A matching source with an invalid key must fail before path use.
        state = json.loads(manifest.read_text())
        state["source"] = genome.genome_fasta_urls[0]
        manifest.write_text(json.dumps(state))
        release()


def test_prune_ignores_unowned_entries_and_external_symlinks(shared_cache, tmp_path):
    genome = release()
    genome.download_genome_fasta()
    shutil.rmtree(genome.download_cache.cache_directory_path)
    objects = genome._genome_fasta.directory.parent
    foreign = objects / ("a" * 64)
    foreign.mkdir()
    (foreign / "precious.fa").write_bytes(DNA)
    external = tmp_path / "external"
    external.mkdir()
    (external / "precious.fa").write_bytes(DNA)
    (objects / ("b" * 64)).symlink_to(external, target_is_directory=True)
    pruned = prune_genome_fastas()
    assert len(pruned) == 1
    assert (foreign / "precious.fa").read_bytes() == DNA
    assert (external / "precious.fa").read_bytes() == DNA


def test_prune_ignores_symlinked_semantic_directories(shared_cache, tmp_path):
    genome = release()
    genome.download_genome_fasta()
    shutil.rmtree(genome.download_cache.cache_directory_path)
    species_directory = cache.dna_cache_root() / "homo_sapiens"
    external = tmp_path / "moved_species_cache"
    species_directory.rename(external)
    species_directory.symlink_to(external, target_is_directory=True)
    assert prune_genome_fastas() == []
    assert Path(genome.genome_fasta_path).read_bytes() == DNA


def test_misplaced_object_descriptor_does_not_authorize_deletion(
    shared_cache, tmp_path
):
    genome = release()
    genome.download_genome_fasta()
    shutil.rmtree(genome.download_cache.cache_directory_path)
    foreign = genome._genome_fasta.directory.parent / "my-data"
    shutil.copytree(genome._genome_fasta.directory, foreign)
    pruned = prune_genome_fastas()
    assert [path for path, _ in pruned] == [str(genome._genome_fasta.directory)]
    assert (foreign / "sequence.fa").read_bytes() == DNA


def test_short_key_collision_does_not_reuse_or_overwrite_data(
    shared_cache, monkeypatch
):
    original_key = cache._identity_key
    monkeypatch.setattr(
        cache, "_identity_key", lambda identity: "a" * 16 + original_key(identity)[16:]
    )
    first = release(81)
    first.download_genome_fasta()
    original_metadata = cache.urlopen

    def changed(request, **kwargs):
        response = original_metadata(request, **kwargs)
        if isinstance(request, str) and request.endswith("CHECKSUMS"):
            response = io.BytesIO(response.read().replace(b"17918", b"17919"))
        return response

    monkeypatch.setattr(cache, "urlopen", changed)
    second = release(82)
    with pytest.raises(ValueError, match="Conflicting shared genome FASTA identity"):
        second.download_genome_fasta()
    assert len(shared_cache) == 1
    assert first.sequence("MT", 1, 4) == "GCTA"
    assert not second._genome_fasta.reference_path.exists()
    first.close()


@pytest.mark.parametrize(
    "field,value",
    [
        ("species", "../outside"),
        ("assembly", "../../GCA_000001405.18"),
        ("provider", "../outside"),
        ("filename", "Homo_sapiens.../outside.dna.toplevel.fa.gz"),
    ],
)
def test_invalid_semantic_path_components_are_rejected(
    shared_cache, monkeypatch, field, value
):
    original_identity = cache._remote_identity

    def invalid(source):
        identity = original_identity(source)
        identity[field] = value
        return identity

    monkeypatch.setattr(cache, "_remote_identity", invalid)
    genome = release()
    with pytest.raises(ValueError, match="Invalid .+ in shared genome FASTA identity"):
        genome.download_genome_fasta()
    assert shared_cache == []
    assert not genome._genome_fasta.reference_path.exists()


def test_failed_shared_download_leaves_only_prunable_owned_object(
    shared_cache, monkeypatch, tmp_path
):
    serve_downloads(monkeypatch, tmp_path, lambda url: b"bad FASTA")
    genome = release()
    with pytest.raises(ValueError, match="FASTA header|size mismatch"):
        genome.download_genome_fasta()
    assert not genome._genome_fasta.reference_path.exists()
    assert not genome._genome_fasta.manifest_path.exists()
    assert len(prune_genome_fastas()) == 1


def test_download_blocks_neither_other_objects_nor_prune(
    shared_cache, monkeypatch, tmp_path
):
    installed = release(81, genome_fasta_mask="soft")
    installed.download_genome_fasta()
    installed.close()
    entered = threading.Event()
    finish_download = threading.Event()
    errors = []

    def slow_payload(url):
        entered.set()
        assert finish_download.wait(10)
        return gzip.compress(DNA)

    serve_downloads(monkeypatch, tmp_path, slow_payload)
    downloading = release(81)

    def install():
        try:
            downloading.download_genome_fasta()
        except Exception as error:
            errors.append(error)

    def meanwhile():
        # A long download must not hold any lock other work needs.
        try:
            other = release(81, genome_fasta_mask="soft")
            other.index_genome_fasta(overwrite=True)
            assert other.sequence("MT", 1, 4) == "GCTA"
            other.close()
            # The unregistered object being downloaded is busy, not orphaned.
            assert prune_genome_fastas() == []
        except BaseException as error:
            errors.append(error)

    worker = threading.Thread(target=install)
    worker.start()
    try:
        assert entered.wait(10)
        checker = threading.Thread(target=meanwhile)
        checker.start()
        checker.join(10)
        assert not checker.is_alive(), "work waited for an unrelated download"
    finally:
        finish_download.set()
        worker.join(10)
    assert not worker.is_alive() and errors == []
    assert downloading.sequence("MT", 1, 4) == "GCTA"
    assert prune_genome_fastas() == []
    downloading.close()


def test_shared_reader_uses_same_downstream_fasta_protocol(shared_cache):
    genome = release()
    genome.download_genome_fasta()
    with genome:
        # Varcode's reference_range accesses .fasta[contig][start-1:end].seq.
        assert genome.fasta["1"][2:10].seq.upper() == "GTNNTTAA"
        restored = Genome.from_dict(Genome("x", "y").to_dict())
        assert restored.fasta is None


def test_unreadable_release_directory_aborts_pruning(shared_cache, monkeypatch):
    genome = release()
    genome.download_genome_fasta()
    original_iterdir = Path.iterdir
    annotation_dir = Path(genome.download_cache.cache_directory_path)

    def unreadable(path):
        if path == annotation_dir / "genome_fasta_refs":
            raise PermissionError("unreadable references")
        return original_iterdir(path)

    monkeypatch.setattr(Path, "iterdir", unreadable)
    with pytest.raises(ValueError, match="unable to inspect"):
        prune_genome_fastas()
    assert Path(genome.genome_fasta_path).exists()


def test_shared_index_rebuild_after_another_instance_overwrites(shared_cache):
    first = release()
    first.download_genome_fasta()
    assert first.sequence("MT", 1, 4) == "GCTA"
    old_reader = first.fasta
    second = release()
    second.download_genome_fasta(overwrite=True)
    assert first.sequence("MT", 1, 4) == "GCTA"
    assert first.fasta is not old_reader
    # Atomic replacement keeps a reader held elsewhere consistent.
    assert old_reader["MT"][0:4].seq == "GCTA"
    first.close()
    second.close()


def test_download_must_match_recorded_identity_size(shared_cache, monkeypatch):
    original = cache.urlopen

    def wrong_size(request, **kwargs):
        response = original(request, **kwargs)
        if not isinstance(request, str):
            response.headers["Content-Length"] = str(
                int(response.headers["Content-Length"]) + 1
            )
        return response

    monkeypatch.setattr(cache, "urlopen", wrong_size)
    genome = release()
    with pytest.raises(ValueError, match="size mismatch"):
        genome.download_genome_fasta()
    assert genome.genome_fasta_path is None
    assert not genome._genome_fasta.reference_path.exists()


def installed(number=81, **kwargs):
    genome = release(number, **kwargs)
    genome.download_genome_fasta()
    genome.index_genome_fasta()
    genome.close()
    return genome


def set_tree_writable(root, writable):
    for current, _, files in os.walk(root):
        os.chmod(current, 0o755 if writable else 0o555)
        for name in files:
            os.chmod(os.path.join(current, name), 0o644 if writable else 0o444)


def test_ensembl_genomes_identity_shares_dna_across_releases(shared_cache):
    genomes = [
        EnsemblRelease(number, species="arabidopsis_thaliana", download_genome_fasta=True)
        for number in (57, 58)
    ]
    for genome in genomes:
        genome.download_genome_fasta()
    assert len(shared_cache) == 1
    assert genomes[0].genome_fasta_path == genomes[1].genome_fasta_path
    assert "TAIR10-GCA_000001735.1" in Path(genomes[0].genome_fasta_path).parts
    assert genomes[1].sequence("MT", 1, 4) == "GCTA"
    genomes[1].close()


def test_published_files_follow_the_umask(shared_cache):
    previous = os.umask(0o002)
    try:
        genome = installed()
    finally:
        os.umask(previous)
    dna = genome._genome_fasta
    for path in (
        dna.materialized_path,
        dna.index_path,
        dna.index_state_path,
        dna.directory / "object.json",
        dna.reference_path,
        dna.manifest_path,
    ):
        assert stat.S_IMODE(path.stat().st_mode) == 0o664, path


def test_installed_cache_serves_readers_without_write_access(shared_cache, tmp_path):
    installed()
    cache_root = tmp_path / "pyensembl"
    set_tree_writable(cache_root, False)
    try:
        genome = release(81)
        assert genome.sequence("MT", 1, 4) == "GCTA"
        assert genome.fasta["1"][2:10].seq == "gtNNTTaa"
        # Repeating completed installation steps writes nothing.
        genome.download_genome_fasta()
        genome.index_genome_fasta()
        genome.close()
    finally:
        set_tree_writable(cache_root, True)


def test_warm_lookups_only_stat_files(shared_cache, monkeypatch):
    genome = release(81)
    genome.download_genome_fasta()
    assert genome.sequence("MT", 1, 4) == "GCTA"

    def unexpected(*args, **kwargs):
        pytest.fail("A warm lookup re-resolved files, read JSON, or locked")

    monkeypatch.setattr(GenomeFasta, "prepare", unexpected)
    monkeypatch.setattr(cache.SharedGenomeFasta, "prepare", unexpected)
    monkeypatch.setattr("pyensembl.genome_fasta._read_json", unexpected)
    monkeypatch.setattr(cache, "_read_json", unexpected)
    monkeypatch.setattr(cache, "FileLock", unexpected)
    for _ in range(3):
        assert genome.sequence("1", 3, 10) == "GTNNTTAA"
        assert genome.fasta["MT"][0:4].seq == "GCTA"
    genome.close()


def test_release_private_copy_never_stands_in_for_shared_dna(shared_cache):
    shared = release(81)
    private = Genome(
        "GRCh38",
        "ensembl",
        annotation_version=81,
        genome_fasta_path_or_url=shared.genome_fasta_urls[0],
    )
    assert private.download_cache.cache_directory_path == (
        shared.download_cache.cache_directory_path
    )
    private.download_genome_fasta()
    assert type(private._genome_fasta) is GenomeFasta
    assert shared.fasta is None and shared.genome_fasta_path is None
    with pytest.raises(MissingGenomeFastaError):
        shared.sequence("MT", 1, 4)
    assert not shared._genome_fasta.reference_path.exists()
    shared.download_genome_fasta()
    assert release(81).sequence("MT", 1, 4) == "GCTA"
    assert prune_genome_fastas() == []


def test_symlinked_dna_cache_root_is_supported(
    shared_cache, tmp_path, monkeypatch, capsys
):
    moved = tmp_path / "bigger_disk" / "dna_cache"
    moved.mkdir(parents=True)
    root = cache.dna_cache_root()
    root.parent.mkdir(parents=True)
    root.symlink_to(moved, target_is_directory=True)
    genome = installed()
    assert Path(genome.genome_fasta_path).resolve().is_relative_to(moved.resolve())
    assert release(81).sequence("MT", 1, 4) == "GCTA"
    run_cli(monkeypatch, "list")
    assert "Genome FASTA: downloaded, indexed" in capsys.readouterr().out
    assert prune_genome_fastas() == []
    run_cli(monkeypatch, "delete-all-files", "--release", "81")
    assert len(prune_genome_fastas()) == 1
    assert not list(moved.rglob("sequence.fa"))


def test_list_reports_bad_reference_without_hiding_other_releases(
    shared_cache, monkeypatch, capsys
):
    first, _ = installed(81), installed(82)
    reference = first._genome_fasta.reference_path
    state = json.loads(reference.read_text())
    state["shared_key"] = "a" * 64
    reference.write_text(json.dumps(state))
    run_cli(monkeypatch, "list")
    lines = capsys.readouterr().out.splitlines()
    statuses = {
        line.split("release=")[1].split(",")[0]: lines[number + 1]
        for number, line in enumerate(lines)
        if "EnsemblRelease(release=" in line
    }
    assert "invalid reference" in statuses["81"]
    assert "downloaded, indexed" in statuses["82"]


def test_unnested_cache_layout_keeps_dna_release_private(
    shared_cache, tmp_path, monkeypatch
):
    # Windows appdirs defaults do not place release caches beside dna_cache,
    # where pruning looks for references.
    elsewhere = tmp_path / "elsewhere" / "pyensembl" / "dna_cache"
    monkeypatch.setattr(cache, "dna_cache_root", lambda: elsewhere)
    genome = release(81)
    assert type(genome._genome_fasta) is GenomeFasta
    genome.download_genome_fasta()
    assert genome.sequence("MT", 1, 4) == "GCTA"
    assert not elsewhere.exists()
    genome.close()
