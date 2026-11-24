"""Cross-module public API contracts and private-helper boundary tests (#519)."""

import ast
from pathlib import Path

import pyensembl
from pyensembl import DownloadCache, EnsemblRelease, Genome
from pyensembl.assembly import assembly_matches
from pyensembl.common import (
    atomic_output,
    publish_file,
    read_json,
    remove_if_exists,
    staging_path,
    write_json,
)
from pyensembl.dated_releases import cached_annotation_dates
from pyensembl.download_cache import default_cache_directory, same_file_bytes
from pyensembl.species import Species, find_species_by_name
from pyensembl.versioned_ids import split_versioned_id

PACKAGE_DIR = Path(pyensembl.__file__).resolve().parent


def test_no_cross_module_private_imports():
    """Modules in pyensembl/ must not import single-underscore names from siblings."""
    violations = []
    for path in sorted(PACKAGE_DIR.rglob("*.py")):
        rel = path.relative_to(PACKAGE_DIR)
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.ImportFrom):
                continue
            for alias in node.names:
                name = alias.name
                if name.startswith("_") and not (name.startswith("__") and name.endswith("__")):
                    module = "." * node.level + (node.module or "")
                    violations.append(f"{rel}:{node.lineno}: from {module} import {name}")
    assert violations == []


def test_cli_package_defines_no_private_helpers():
    """Every function and class in pyensembl/cli/ must be public with a docstring (#504)."""
    cli_dir = PACKAGE_DIR / "cli"
    violations = []
    for path in sorted(cli_dir.glob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                if node.name.startswith("_"):
                    violations.append(f"{path.name}:{node.lineno}: private {node.name}")
                if not ast.get_docstring(node):
                    violations.append(f"{path.name}:{node.lineno}: missing docstring on {node.name}")
    assert violations == []


def test_no_cross_module_species_registry_or_genome_private_access():
    """Other modules must use public Species/Genome/DownloadCache APIs instead of private attrs."""
    forbidden_species_attrs = {
        "_latin_names_to_species",
        "_common_names_to_species",
        "_reference_names_to_species",
    }
    forbidden_cross_module_attrs = {
        "_cache_path",
        "_annotation_sources",
        "_annotation_status",
        "_uses_annotation_checksums",
        "_source_translation_table",
        "_reference_metadata",
        "_gtf_path_or_url",
    }
    violations = []
    for path in sorted(PACKAGE_DIR.rglob("*.py")):
        rel = path.relative_to(PACKAGE_DIR)
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Attribute):
                continue
            if isinstance(node.value, ast.Name) and node.value.id in ("self", "cls"):
                continue
            if path.name != "species.py" and node.attr in forbidden_species_attrs:
                violations.append(f"{rel}:{node.lineno}: {ast.unparse(node)}")
            if path.name not in ("genome.py", "ensembl_release.py", "download_cache.py"):
                if node.attr in forbidden_cross_module_attrs:
                    violations.append(f"{rel}:{node.lineno}: {ast.unparse(node)}")
                if node.attr == "_resolved":
                    violations.append(f"{rel}:{node.lineno}: {ast.unparse(node)}")
    assert violations == []


def test_species_public_registry_accessors():
    all_species = Species.all()
    assert all(isinstance(s, Species) for s in all_species)
    assert len(all_species) == len(Species.all_registered_latin_names())

    human = Species.get_by_latin_name("homo_sapiens")
    assert human is find_species_by_name("human")
    assert Species.get_by_latin_name("not_a_real_species") is None

    refs = Species.all_reference_names()
    assert "GRCh38" in refs
    assert "GRCh37" in refs
    assert Species.get_by_reference_name("GRCh38") is human
    assert Species.get_by_reference_name("not_a_real_reference") is None

    assert human.newest_release == pyensembl.MAX_ENSEMBL_RELEASE
    arabidopsis = find_species_by_name("arabidopsis")
    assert arabidopsis.newest_release == max(
        end for _, end in arabidopsis.reference_assemblies.values()
    )


def test_common_file_and_json_helpers(tmp_path):
    target = tmp_path / "data.json"
    assert read_json(target) is None
    write_json(target, {"a": 1, "b": [2, 3]})
    assert read_json(target) == {"a": 1, "b": [2, 3]}

    staged = staging_path(tmp_path, "out.bin")
    assert staged.parent == tmp_path
    assert staged.name.startswith(".out.bin.")
    staged.write_bytes(b"hello")
    final = tmp_path / "out.bin"
    publish_file(staged, final)
    assert final.read_bytes() == b"hello"
    assert not staged.exists()

    with atomic_output(tmp_path / "atomic.txt", "w") as handle:
        handle.write("ok")
    assert (tmp_path / "atomic.txt").read_text() == "ok"

    remove_if_exists(final)
    assert not final.exists()
    remove_if_exists(final)  # No error when already missing.


def test_download_cache_and_versioned_id_public_helpers(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    directory = default_cache_directory("GRCh38", "ensembl", 116)
    assert directory.endswith("GRCh38/ensembl116")

    cache = DownloadCache("GRCh38", "ensembl", 116, cache_directory_path=str(tmp_path))
    url = "https://example.org/data/sample.gtf.gz"
    modern = cache.cache_path(url)
    legacy = cache.cache_path(url, legacy=True)
    assert modern != legacy
    assert legacy.endswith("sample.gtf.gz")

    first = tmp_path / "a.txt"
    second = tmp_path / "b.txt"
    first.write_bytes(b"same bytes")
    second.write_bytes(b"same bytes")
    assert same_file_bytes(str(first), str(second))
    second.write_bytes(b"different")
    assert not same_file_bytes(str(first), str(second))

    assert split_versioned_id("ENST00000269305.8") == ("ENST00000269305", 8)
    assert split_versioned_id("ENST00000269305") == ("ENST00000269305", None)
    assert split_versioned_id("AT1G01010.1") == ("AT1G01010.1", None)


def test_assembly_and_dated_release_public_helpers(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    record = {
        "assembly": "GRCh38",
        "assembly_name": "GRCh38.p14",
        "assembly_accession": "GCA_000001405.29",
        "aliases": ["hg38"],
    }
    assert assembly_matches(record, None)
    assert assembly_matches(record, "grch38")
    assert assembly_matches(record, "GCA_000001405.29")
    assert assembly_matches(record, "hg38")
    assert not assembly_matches(record, "GRCh37")

    assert cached_annotation_dates("GCA_000001405.29", "ensembl") is None


def test_genome_and_ensembl_release_public_helpers(tmp_path):
    genome = Genome(
        "GRCh38", "test", annotation_version="1",
        gtf_path_or_url=str(tmp_path / "missing.gtf"),
        cache_directory_path=str(tmp_path / "cache"),
    )
    assert genome.gtf_path_or_url == str(tmp_path / "missing.gtf")
    assert genome.annotation_sources() == {"gtf": str(tmp_path / "missing.gtf")}
    assert genome.annotation_status() == "missing"
    assert genome.genome_fasta_status() is None
    assert not genome.uses_annotation_checksums()
    assert genome.pinned_source_checksums() == {}

    release = EnsemblRelease("2026_04", cache_directory_path=str(tmp_path / "dated"))
    assert release.uses_annotation_checksums()
    assert release.pinned_source_checksums() == {}
    release.pin_source_checksum(
        "https://example.org/genes.gtf.gz",
        {"source": "https://example.org/genes.gtf.gz", "md5": "0" * 32},
    )
    assert release.pinned_source_checksums() == {
        "https://example.org/genes.gtf.gz": {
            "source": "https://example.org/genes.gtf.gz",
            "md5": "0" * 32,
        }
    }


def test_annotation_integrity_acquire_and_verify_use_public_genome_apis(tmp_path, monkeypatch):
    """Exercise annotation_integrity.acquire/verify_sources offline without sockets."""
    import gzip
    import hashlib
    import datacache
    from pyensembl import annotation_integrity

    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    release = EnsemblRelease("2026_04", cache_directory_path=str(tmp_path / "dated"))
    payload = gzip.compress(b"1\ttest\tgene\t1\t10\t.\t+\t.\tgene_id \"ENSG1\";\n")
    md5 = hashlib.md5(payload, usedforsecurity=False).hexdigest()
    manifest = "".join(
        f"{md5}  {url.rsplit('/', 1)[1]}\n"
        for url in release.annotation_sources().values()
    ).encode("utf-8")

    monkeypatch.setattr(datacache, "fetch_bytes", lambda url, **kw: manifest)

    def fake_fetch(self, url, destination, overwrite, show_progress=False, validator=None):
        Path(destination).parent.mkdir(parents=True, exist_ok=True)
        Path(destination).write_bytes(payload)
        if validator is not None:
            validator(destination)
        return destination

    monkeypatch.setattr(DownloadCache, "_fetch", fake_fetch)

    path = annotation_integrity.acquire(release, "gtf", release.gtf_url, False, False)
    assert Path(path).is_file()
    assert release.gtf_url in release.pinned_source_checksums()
    assert release.pinned_source_checksums()[release.gtf_url]["md5"] == md5

    reports = {r["role"]: r for r in release.verify_sources(check_remote=True)}
    assert reports["gtf"]["status"] == "verified"
    assert reports["gtf"]["upstream_status"] == "unchanged"

