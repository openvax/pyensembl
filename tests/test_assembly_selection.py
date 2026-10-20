"""Assembly selection against published metadata and pre-change saved objects."""

import base64
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import pickle
import shlex
import tarfile
from unittest.mock import patch

import pytest

from pyensembl import EnsemblRelease, Genome, shell
from pyensembl import assembly as assemblies
from pyensembl.species import Species


DATA = Path(__file__).parent / "data" / "assembly_metadata"


@pytest.fixture(autouse=True)
def isolated_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path / "cache"))


@pytest.mark.parametrize("species, release", [
    ("homo_sapiens", 75), ("homo_sapiens", 93), ("mus_musculus", 81),
    ("canis_familiaris", 104), ("canis_familiaris", 105),
    ("canis_lupus_familiarisbasenji", 116), ("arabidopsis_thaliana", 63),
])
def test_numbered_metadata_matches_captured_upstream_header(species, release):
    capture = json.loads((DATA / ("%s-%s.json" % (species, release))).read_text())
    fields = dict(line[2:].split(" ", 1) for line in capture["headers"])
    genome = EnsemblRelease(release, species=species)
    assert genome.assembly == fields["genome-version"]
    assert genome.assembly_name == fields["genome-build"]
    assert genome.assembly_accession == fields["genome-build-accession"].rsplit(":", 1)[-1]


def test_all_registered_releases_construct_offline_with_metadata():
    # Audit coverage independently through the public species registry. Missing
    # upstream accessions are documented explicitly, not silently fabricated.
    with patch("urllib.request.urlopen", side_effect=AssertionError("network forbidden")):
        for species, release in Species.all_species_release_pairs():
            genome = EnsemblRelease(release, species=species)
            assert genome.assembly
            assert genome.assembly_name
            assert hasattr(genome, "assembly_accession")
            if genome.assembly_accession is None:
                assert any(
                    gap["species"] == species and gap["release"] == release
                    for gap in assemblies.catalogue()["metadata_gaps"]
                )


@pytest.mark.parametrize("selector", [None, "GRCh37", "grch37", "GRCh37.p13", "GCA_000001405.14"])
def test_historical_dated_selection_infers_full_identity(selector):
    options = {} if selector is None else dict(assembly=selector)
    genome = EnsemblRelease("2013_09", **options)
    assert (genome.assembly, genome.assembly_name, genome.assembly_accession) == (
        "GRCh37", "GRCh37.p13", "GCA_000001405.14",
    )
    assert genome.species.latin_name == "homo_sapiens"
    assert "/GCA/000/001/405/14/ensembl/2013_09/" in genome.gtf_url


@pytest.mark.parametrize("release, options", [
    (75, dict(assembly="GRCh38")),
    (93, dict(assembly="GRCh38.p14")),
    (75, dict(assembly_accession="GCA_000001405.29")),
    ("2013_09", dict(assembly="GRCh37", species="mouse")),
    ("2013_09", dict(assembly="GRCh37", assembly_accession="GCA_000001405.29")),
    ("2019_08", dict(assembly="Basenji_breed-1.1", species="dog_great_dane")),
])
def test_conflicting_selectors_fail_without_cache_writes(release, options, tmp_path):
    with pytest.raises(ValueError):
        EnsemblRelease(release, **options)
    assert not (tmp_path / "cache").exists()


def test_numbered_named_assertion_keeps_existing_cache():
    plain = EnsemblRelease(75)
    named = EnsemblRelease(75, assembly="GRCh37.p13")
    assert named.download_cache.cache_directory_path == plain.download_cache.cache_directory_path
    assert named.reference_name == "GRCh37"


def test_species_default_and_accession_compatibility_layouts():
    plain = EnsemblRelease("2026_04")
    named = EnsemblRelease("2026_04", assembly="GRCh38.p14")
    explicit = EnsemblRelease("2026_04", assembly_accession="GCA_000001405.29")
    assert plain.reference_name == named.reference_name == "GRCh38"
    assert plain.include_alt is named.include_alt is True
    assert explicit.reference_name == "GCA_000001405.29"
    assert explicit.include_alt is False
    assert plain.download_cache.cache_directory_path == named.download_cache.cache_directory_path
    assert plain.download_cache.cache_directory_path != explicit.download_cache.cache_directory_path
    for selector in ["GCA_000001405.29", "gca_000001405.29", " GCA_000001405.29 "]:
        selected = EnsemblRelease("2026_04", assembly=selector)
        assert selected.gtf_url == explicit.gtf_url
        assert selected.download_cache.cache_directory_path == explicit.download_cache.cache_directory_path


def test_existing_dated_directory_and_mirror_keep_original_selection(tmp_path):
    old = EnsemblRelease.from_dict(dict(release="2013_09", species="human"))
    assert old.assembly == "GRCh38"
    directory = Path(old.download_cache.cache_directory_path)
    directory.mkdir(parents=True)
    marker = directory / "existing"
    marker.write_bytes(b"leave alone")
    before = marker.stat()
    reopened = EnsemblRelease("2013_09")
    assert reopened.assembly == "GRCh38"
    assert reopened.download_cache.cache_directory_path == str(directory)
    assert marker.read_bytes() == b"leave alone"
    assert marker.stat().st_mtime_ns == before.st_mtime_ns
    assert marker.stat().st_ino == before.st_ino
    assert EnsemblRelease("2013_09", server="https://mirror.example").assembly == "GRCh38"


def test_dogs_keep_reference_breeds_and_meaningful_versions():
    old = EnsemblRelease(104, species="dog")
    current = EnsemblRelease("2020_10", species="dog")
    basenji = EnsemblRelease("2019_08", assembly="Basenji_breed-1.1")
    assert old.assembly == old.assembly_name == "CanFam3.1"
    assert current.assembly_name == "ROS_Cfam_1.0"
    assert basenji.species.latin_name == "canis_lupus_familiarisbasenji"
    assert basenji.assembly_name == "Basenji_breed-1.1"
    assert basenji.assembly_accession == "GCA_004886185.1"
    assert len({old.assembly_accession, current.assembly_accession, basenji.assembly_accession}) == 3
    assert EnsemblRelease("2019_08", species="dog", assembly=basenji.assembly_name).assembly_accession == basenji.assembly_accession
    with pytest.raises(ValueError, match="Ambiguous"):
        EnsemblRelease("2019_08", species="dog")


def test_fly_exact_metadata_does_not_rename_legacy_reference():
    genome = EnsemblRelease("2022_07", species="fly")
    assert genome.reference_name == "BDGP6.54"
    assert genome.assembly_name == "BDGP6.46"
    assert "BDGP6.54" in genome.download_cache.cache_directory_path


def test_unknown_accession_does_not_become_human():
    genome = EnsemblRelease("2023_03", assembly="GCA_999999999.1")
    assert genome.species is None
    assert genome.assembly is None
    assert genome.assembly_name is None
    assert genome.assembly_accession == "GCA_999999999.1"
    assert "GCA_999999999.1" in genome.install_string()
    assert EnsemblRelease.from_json(genome.to_json()).species is None


@pytest.mark.parametrize("serializer", ["json", "pickle"])
def test_published_2283_release_objects_still_load(serializer):
    for state in json.loads((DATA / "pyensembl-2.28.3.json").read_text()):
        if json.loads(state["json"])["__class__"]["__name__"] != "EnsemblRelease":
            continue  # The explicitly retired class no longer has a loader.
        if serializer == "json":
            genome = Genome.from_json(state["json"])
        else:
            genome = pickle.loads(base64.b64decode(state["pickle"]))
        assert genome.reference_name == state["reference_name"]
        assert genome.assembly_name
        assert genome.assembly_accession
        restored = type(genome).from_json(genome.to_json())
        assert restored == genome
        assert restored.download_cache.cache_directory_path == genome.download_cache.cache_directory_path


def test_serialized_selection_survives_catalogue_changes(monkeypatch):
    genome = EnsemblRelease("2026_04", genome_fasta=True)
    encoded = genome.to_json()
    payload = pickle.dumps(genome)
    before = (genome.assembly, genome.assembly_name, genome.assembly_accession,
              genome.reference_name, genome.gtf_url, genome.download_cache.cache_directory_path)
    monkeypatch.setattr(assemblies, "catalogue", lambda: dict(numbered={}, dated=[]))
    for restored in [EnsemblRelease.from_json(encoded), pickle.loads(payload)]:
        assert (restored.assembly, restored.assembly_name, restored.assembly_accession,
                restored.reference_name, restored.gtf_url,
                restored.download_cache.cache_directory_path) == before
        assert restored.requires_genome_fasta


def test_cache_keys_separate_coverage_and_provider():
    primary = EnsemblRelease.cached("2026_04", assembly_accession="GCA_000001405.29")
    alt = EnsemblRelease.cached("2026_04", assembly_accession="GCA_000001405.29", include_alt=True)
    other = EnsemblRelease.cached("2026_04", assembly_accession="GCA_000001405.29", provider="community")
    assert primary is EnsemblRelease.cached("2026_04", assembly_accession="GCA_000001405.29")
    assert len({primary, alt, other}) == 3
    named_other = EnsemblRelease(
        "2026_04", assembly="GRCh38.p14",
        assembly_accession="GCA_000001405.29", provider="community",
    )
    assert named_other.assembly_name == other.assembly_name == "GRCh38.p14"
    assert named_other.gtf_url == other.gtf_url


@pytest.mark.parametrize("release, options", [
    (75, {}),
    ("2013_09", dict(assembly_accession="GCA_000001405.14")),
    ("2023_03", dict(assembly_accession="GCA_999999999.1")),
])
def test_install_hints_preserve_selected_mirror_and_custom_cache(release, options, tmp_path):
    genome = EnsemblRelease(
        release, **options, server="https://mirror.example/root/",
        reference_name="my label", cache_directory_path=str(tmp_path / "custom cache"),
        genome_fasta=True,
    )
    args = shell.parser.parse_args(shlex.split(genome.install_string())[1:])
    reopened, = shell.collect_selected_genomes(args)
    assert reopened.gtf_url == genome.gtf_url
    assert reopened.genome_fasta_urls == genome.genome_fasta_urls
    assert reopened.reference_name == genome.reference_name
    assert reopened.download_cache.cache_directory_path == genome.download_cache.cache_directory_path
    assert reopened.assembly_accession == genome.assembly_accession
    assert "server='https://mirror.example/root/'" in genome._genome_fasta_setup_hint()
    if release == 75:
        # Exact URL construction from the published numbered-release helper:
        # the pre-existing double slash must not become a new source identity.
        assert "/root//pub/" in genome.gtf_url


def test_ambiguous_nondefault_assemblies_are_not_guessed(monkeypatch):
    data = deepcopy(assemblies.catalogue())
    old = next(r for r in data["dated"] if r["assembly_accession"] == "GCA_000001405.14")
    second = dict(old, assembly_name="another assembly", assembly_accession="GCA_999999999.1")
    data["dated"].append(second)
    monkeypatch.setattr(assemblies, "catalogue", lambda: data)
    with pytest.raises(ValueError, match="Ambiguous"):
        EnsemblRelease("2013_09")
    assert EnsemblRelease("2013_09", assembly="GRCh37.p13").assembly_accession == old["assembly_accession"]


def test_cli_uses_shared_assembly_resolution():
    args = shell.parser.parse_args([
        "inspect", "--release", "2013_09", "--assembly", "GRCh37.p13",
    ])
    genome, = shell.collect_selected_genomes(args)
    assert genome.assembly_accession == "GCA_000001405.14"
    args = shell.parser.parse_args([
        "inspect", "--release", "2013_09", "--reference-name", "GRCh37",
    ])
    assert shell.collect_selected_genomes(args)[0].assembly_name == "GRCh37.p13"


def test_cli_rejects_ambiguous_cache_and_species_combinations():
    for options in [
        ["--species", "human", "mouse", "--assembly", "GRCh37"],
        ["--cache-directory", "/unused", "--release", "75", "76"],
    ]:
        args = shell.parser.parse_args(["install", "--release", "75", *options])
        with pytest.raises(ValueError):
            shell.collect_selected_genomes(args)


@pytest.mark.parametrize("read_only", [False, True])
def test_published_installed_caches_reopen_without_mutation(tmp_path, monkeypatch, read_only):
    from pyensembl import Database

    root = tmp_path / "published"
    with tarfile.open(DATA / "legacy-cache-2.28.3.tar.gz") as archive:
        archive.extractall(root, filter="data")
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(root / "cache"))

    def forbidden(*args, **kwargs):
        pytest.fail("Existing installed cache attempted network access or parsing")

    monkeypatch.setattr("requests.sessions.Session.request", forbidden)
    monkeypatch.setattr("urllib.request.urlopen", forbidden)
    monkeypatch.setattr("pyensembl.sequence_data.parse_fasta_dictionary", forbidden)
    monkeypatch.setattr(Database, "_load_gtf_as_dataframe", forbidden)
    paths = list(root.rglob("*"))
    modes = {path: path.stat().st_mode for path in paths}
    if read_only:
        for path in paths:
            path.chmod(0o555 if path.is_dir() else 0o444)

    def snapshot():
        result = {}
        for path in root.rglob("*"):
            if path.is_file():
                st = path.stat()
                result[str(path.relative_to(root))] = (
                    hashlib.sha256(path.read_bytes()).hexdigest(), st.st_mode,
                    st.st_ino, st.st_mtime_ns, st.st_ctime_ns,
                )
        return result

    before = snapshot()
    try:
        for state in json.loads((DATA / "legacy-cache-objects.json").read_text()):
            # Only relocate the explicitly chosen custom cache; source URLs
            # remain exactly as saved by the published version.
            encoded = state["json"].replace(
                "/private/tmp/pyensembl-482-legacy/custom-cache", str(root / "custom-cache")
            )
            settings = json.loads(encoded)
            if settings["__class__"]["__name__"] == "EnsemblRelease":
                genome = Genome.from_json(encoded)
            else:
                # Reopen files created by the retired constructor through the
                # supported API; keep the historical fixture itself unchanged.
                settings.pop("__class__")
                genome = EnsemblRelease(settings.pop("annotation_date"), **settings)
            with genome:
                genome.download()
                genome.index()
                assert genome.transcript_by_id("ENST1").sequence == "CCCCGGTT"
                assert genome.protein_sequences.get("ENSP1.1") == "MP"
                assert genome.reference_name == state["reference_name"]
        assert snapshot() == before
    finally:
        if read_only:
            for path, mode in modes.items():
                path.chmod(mode)
