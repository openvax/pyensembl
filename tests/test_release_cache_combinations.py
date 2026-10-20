"""Install real tiny datasets through different EnsemblRelease selections."""

from dataclasses import dataclass
import gzip
import hashlib
from io import BytesIO
import json
from pathlib import Path
import pickle
import shlex

import pytest
import requests

from pyensembl import Database, EnsemblRelease, shell
from pyensembl import ensembl_metadata, genome_fasta_cache
from .test_explicit_ensembl_release import bgzf
from .test_legacy_cache_compatibility import snapshot


HUMAN = "GCA_000001405.29"
GRCH37 = "GCA_000001405.14"
DOG = "GCA_014441545.1"
BASENJI = "GCA_004886185.1"
DANE = "GCA_005444595.1"


@dataclass(frozen=True)
class Selection:
    name: str
    release: object
    marker: str
    accession: str
    variants: tuple
    alt: bool = False


SELECTIONS = [
    Selection("numbered-human", 75, "H75", GRCH37, (
        {}, {"species": "human"}, {"species": "homo_sapiens"},
        {"assembly": "GRCh37"}, {"assembly": "GRCh37.p13"},
        {"assembly": GRCH37}, {"assembly_accession": GRCH37},
        {"assembly": "GRCh37", "assembly_accession": GRCH37, "species": "human"},
    )),
    Selection("dated-human-default", "2026_04", "H26", HUMAN, (
        {}, {"species": "human"}, {"species": "homo_sapiens"},
        {"assembly": "GRCh38"}, {"assembly": "GRCh38.p14"},
        {"include_alt": True}, {"provider": "ensembl"},
    ), alt=True),
    Selection("dated-human-primary", "2026_04", "H26", HUMAN, (
        {"assembly": HUMAN}, {"assembly": " gca_000001405.29 "},
        {"assembly_accession": HUMAN},
        {"assembly_accession": HUMAN, "species": "human"},
        {"assembly_accession": HUMAN, "assembly": "GRCh38.p14", "provider": "ensembl"},
        {"include_alt": False},
    )),
    Selection("dated-human-explicit-alt", "2026_04", "H26", HUMAN, (
        {"assembly": HUMAN, "include_alt": True},
        {"assembly_accession": HUMAN, "include_alt": True},
        {"assembly_accession": HUMAN, "include_alt": True, "species": "human"},
    ), alt=True),
    Selection("historical-human", "2013_09", "H37", GRCH37, (
        {}, {"species": "human"}, {"assembly": "GRCh37"},
        {"assembly": "GRCh37.p13"}, {"assembly": GRCH37},
        {"assembly_accession": GRCH37},
    )),
    Selection("dated-dog-default", "2020_10", "DOG", DOG, (
        {"species": "dog"}, {"species": "canis_familiaris"},
        {"assembly": "ROS_Cfam_1.0"}, {"species": "dog", "include_alt": False},
    )),
    Selection("numbered-dog", 104, "DOG104", "GCA_000002285.2", (
        {"species": "dog"}, {"species": "canis_familiaris"},
        {"assembly": "CanFam3.1"}, {"assembly_accession": "GCA_000002285.2"},
    )),
    Selection("dated-dog-explicit", "2020_10", "DOG", DOG, (
        {"assembly": DOG}, {"assembly_accession": DOG},
        {"assembly_accession": DOG, "species": "dog"},
        {"assembly_accession": DOG, "species": "canis_familiaris", "provider": "ensembl"},
    )),
    Selection("basenji-default", "2019_08", "BASENJI", BASENJI, (
        {"species": "dog_basenji"}, {"assembly": "Basenji_breed-1.1"},
    )),
    Selection("basenji-explicit", "2019_08", "BASENJI", BASENJI, (
        {"assembly": BASENJI}, {"assembly_accession": BASENJI},
        {"assembly_accession": BASENJI, "species": "dog"},
        {"assembly": "Basenji_breed-1.1", "species": "dog"},
    )),
    Selection("great-dane-explicit", "2019_08", "DANE", DANE, (
        {"assembly": DANE}, {"assembly_accession": DANE},
        {"assembly_accession": DANE, "species": "dog"},
        {"assembly": "UMICH_Zoey_3.1", "species": "dog"},
    )),
    Selection("other-provider", "2026_04", "COMMUNITY", HUMAN, (
        {"provider": "community"},
        {"assembly_accession": HUMAN, "provider": "community"},
        {"assembly_accession": HUMAN, "assembly": "GRCh38.p14", "provider": "community"},
    )),
    Selection("other-date", "2025_12", "H25", HUMAN, (
        {"assembly": HUMAN}, {"assembly_accession": HUMAN},
        {"species": "human", "include_alt": False},
    )),
    Selection("custom-label", "2026_04", "H26", HUMAN, (
        {"assembly": HUMAN, "reference_name": "study label"},
        {"assembly_accession": HUMAN, "reference_name": "study label"},
    )),
    Selection("custom-directory", "2026_04", "H26", HUMAN, (
        {"assembly": HUMAN, "cache_directory_path": "{cache}/chosen directory"},
        {"assembly_accession": HUMAN, "cache_directory_path": "{cache}/chosen directory"},
    )),
    Selection("structured-mirror", "2026_04", "MIRROR", HUMAN, (
        {"server": "https://mirror.example/ensembl"},
        {"server": "https://mirror.example/ensembl/", "species": "human"},
        {"server": "https://mirror.example/ensembl", "assembly": "GRCh38.p14"},
    ), alt=True),
    Selection("unknown-accession", "2026_04", "UNKNOWN", "GCF_999999999.1", (
        {"assembly": "GCF_999999999.1"}, {"assembly_accession": "GCF_999999999.1"},
    )),
]


class TinySources:
    """Replace only HTTP transport; use real acquisition, parsers and indexes."""

    def __init__(self, root):
        self.root = root
        self.payloads = {}
        self.downloads = []
        self.offline = False

    def register(self, genome, marker, alt=False):
        gene, transcript, protein = "g_" + marker, "t_" + marker, "p_" + marker
        attrs = (
            'gene_id "%s"; gene_name "%s"; transcript_id "%s"; transcript_version "1";'
            % (gene, marker, transcript)
        )
        rows = (
            '1\tensembl\tgene\t1\t8\t.\t+\t.\tgene_id "%s"; gene_name "%s";\n'
            % (gene, marker)
            + '1\tensembl\ttranscript\t1\t8\t.\t+\t.\t%s\n' % attrs
            + '1\tensembl\texon\t1\t8\t.\t+\t.\t%s exon_id "e"; exon_number "1";\n' % attrs
        )
        if alt:
            rows += 'ALT\tensembl\tgene\t1\t4\t.\t+\t.\tgene_id "alt_%s"; gene_name "%s";\n' % (
                marker, marker,
            )
        digest = hashlib.sha256(marker.encode()).digest()
        cdna = "".join("ACGT"[n % 4] for n in digest[:8])
        peptide = "M" + "ACDEFGHIKLMNPQRSTVWY"[digest[0] % 19]
        expected = dict(gene=gene, transcript=transcript, protein=protein,
                        cdna=cdna, peptide=peptide, alt=alt, marker=marker)

        def publish(url, content):
            payload = bgzf(content) if url.endswith(".bgz") else gzip.compress(content, mtime=0)
            if url in self.payloads:
                assert self.payloads[url] == payload, "one source URL cannot describe two datasets"
            self.payloads[url] = payload

        publish(genome.gtf_url, rows.encode())
        for url in genome.transcript_fasta_urls:
            sequence = ">nc_%s\nGCGC\n" % marker if ".ncrna." in url else (
                ">%s.1\n%s\n" % (transcript, cdna)
            )
            publish(url, sequence.encode())
        for url in genome.protein_fasta_urls:
            publish(url, (">%s.1\n%s\n" % (protein, peptide)).encode())
        for url in genome.genome_fasta_urls:
            # One assembly has the same DNA across dates/coverage/providers;
            # shared-cache identity still distinguishes different providers.
            dna_key = genome.assembly_accession
            if "mirror.example" in url:
                dna_key += "/mirror"
            digest = hashlib.sha256(dna_key.encode()).digest()
            dna = "".join("ACGT"[n % 4] for n in digest[:16])
            if genome.genome_fasta_mask == "soft":
                dna = dna[:4] + dna[4:8].lower() + dna[8:]
            elif genome.genome_fasta_mask == "hard":
                dna = dna[:4] + "NNNN" + dna[8:]
            expected["dna"] = dna
            publish(url, (">1\n%s\n>ALT\nACGT\n" % dna).encode())
        return expected

    def request(self, session, method, url, **kwargs):
        if self.offline:
            pytest.fail("Reopening an installed selection attempted a download: " + url)
        assert method.upper() == "GET" and url in self.payloads, url
        self.downloads.append(url)
        response = requests.Response()
        response.status_code = 200
        response.url = url
        response.headers["Content-Length"] = str(len(self.payloads[url]))
        response.raw = BytesIO(self.payloads[url])
        return response

    def metadata(self, request, **kwargs):
        if self.offline:
            pytest.fail("Reopening an installed selection fetched DNA metadata")
        url = request if isinstance(request, str) else request.full_url
        if url.endswith("/md5sum.txt"):
            directory = url.rsplit("/", 1)[0]
            entries = [
                "%s  %s" % (hashlib.md5(data, usedforsecurity=False).hexdigest(),
                            source.rsplit("/", 1)[1])
                for source, data in self.payloads.items()
                if source.rsplit("/", 1)[0] == directory
            ]
            return BytesIO("\n".join(entries).encode())
        if not isinstance(request, str) and request.get_method() == "HEAD":
            response = BytesIO()
            response.headers = {"Content-Length": str(len(self.payloads[url]))}
            return response
        # Numbered releases use URL-isolated DNA here; dedicated DNA-cache
        # tests exercise their README/CHECKSUMS protocol.
        raise OSError("No numbered-release sharing metadata in this fixture")


@pytest.fixture
def sources(tmp_path, monkeypatch):
    root = tmp_path / "cache"
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(root))
    sources = TinySources(root)
    monkeypatch.setattr(requests.sessions.Session, "request",
                        lambda session, method, url, **kw: sources.request(session, method, url, **kw))
    monkeypatch.setattr(genome_fasta_cache, "urlopen", sources.metadata)
    return sources


def construct(selection, variant, root, factory=EnsemblRelease, **options):
    kwargs = {k: v.replace("{cache}", str(root)) if isinstance(v, str) else v
              for k, v in variant.items()}
    return factory(selection.release, **{"genome_fasta": True, **kwargs, **options})


def query(genome, expected):
    assert genome.installed()
    expected_genes = {expected["gene"]}
    if expected["alt"]:
        expected_genes.add("alt_" + expected["marker"])
    assert set(genome.gene_ids()) == expected_genes
    assert genome.gene_name_of_gene_id(expected["gene"]) == expected["marker"]
    assert genome.transcript_by_id(expected["transcript"]).sequence == expected["cdna"]
    assert genome.protein_sequences.get(expected["protein"] + ".1") == expected["peptide"]
    if genome.requires_genome_fasta:
        assert genome.sequence("1", 1, 16, mask="raw") == expected["dna"]
    return (genome.gtf_path, genome.db.local_db_path,
            tuple(genome.transcript_sequences.fasta_paths),
            tuple(genome.protein_sequences.fasta_paths), genome.genome_fasta_path)


def forbid_parsing(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("Existing cache attempted to regenerate an index")
    monkeypatch.setattr(Database, "_load_gtf_as_dataframe", forbidden)
    monkeypatch.setattr("pyensembl.sequence_data.parse_fasta_dictionary", forbidden)


@pytest.mark.parametrize("selection, initial_options", [
    pytest.param(selection, options, id="%s-%d" % (selection.name, index))
    for selection in SELECTIONS
    for index, options in enumerate(selection.variants)
])
def test_equivalent_arguments_reopen_installed_data_read_only(
    sources, monkeypatch, selection, initial_options,
):
    with construct(selection, initial_options, sources.root) as first:
        assert first.assembly_accession == selection.accession
        expected = sources.register(first, selection.marker, selection.alt)
        assert not sources.root.exists()
        first.download()
        first.index()
        paths = query(first, expected)
        encoded, pickled, command = first.to_json(), pickle.dumps(first), first.install_string()
    source_count = len(sources.downloads)
    sources.offline = True
    forbid_parsing(monkeypatch)
    modes = {path: path.stat().st_mode for path in sources.root.rglob("*")}
    modes[sources.root] = sources.root.stat().st_mode
    for path in modes:
        path.chmod(0o555 if path.is_dir() else 0o444)
    before = snapshot(sources.root)
    try:
        reopeners = [
            lambda variant=variant: construct(selection, variant, sources.root)
            for variant in selection.variants
        ] + [
            lambda: construct(selection, initial_options, sources.root,
                              factory=EnsemblRelease.cached),
            lambda: EnsemblRelease.from_json(encoded),
            lambda: pickle.loads(pickled),
            lambda: shell.collect_selected_genomes(shell.parser.parse_args(shlex.split(command)[1:]))[0],
        ]
        for reopen in reopeners:
            with reopen() as genome:
                genome.download()
                genome.index()
                assert query(genome, expected) == paths
        assert len(sources.downloads) == source_count
        assert snapshot(sources.root) == before
    finally:
        for path, mode in modes.items():
            path.chmod(mode)


def test_distinct_selections_coexist_without_cross_contamination(sources):
    installed = {}
    for selection in SELECTIONS:
        with construct(selection, selection.variants[0], sources.root) as genome:
            expected = sources.register(genome, selection.marker, selection.alt)
            genome.download()
            genome.index()
            paths = query(genome, expected)
            installed[selection.name] = (genome.to_json(), expected, paths)
    assert len({value[2][1] for value in installed.values()}) == len(SELECTIONS)
    # Same verified assembly DNA is shared across date and annotation coverage.
    shared = [installed[name][2][-1] for name in (
        "dated-human-default", "dated-human-primary", "dated-human-explicit-alt", "other-date",
    )]
    assert len(set(shared)) == 1
    assert installed["other-provider"][2][-1] != shared[0]
    sources.offline = True
    before = snapshot(sources.root)
    for encoded, expected, paths in installed.values():
        with EnsemblRelease.from_json(encoded) as genome:
            genome.download()
            genome.index()
            assert query(genome, expected) == paths
    assert snapshot(sources.root) == before


def test_optional_local_and_masked_dna_reuse_the_installed_annotation(sources, tmp_path, monkeypatch):
    selection = next(case for case in SELECTIONS if case.name == "dated-human-primary")
    options = selection.variants[0]
    with construct(selection, options, sources.root, genome_fasta=None) as first:
        expected = sources.register(first, selection.marker)
        first.download()
        first.index()
        annotation_paths = query(first, expected)[:4]
        assert not first.requires_genome_fasta
        assert first.genome_fasta_path is None
    original_files = {
        path: (path.read_bytes(), path.stat()) for path in sources.root.rglob("*")
        if path.is_file()
    }
    local = tmp_path / "local.fa"
    local.write_text(">1\nAACCGGTTAACCGGTT\n")
    second = tmp_path / "second.fa"
    second.write_text(">1\nTTGGCCAATTGGCCAA\n")
    compressed = tmp_path / "local.fa.gz"
    compressed.write_bytes(gzip.compress(b">1\nTTAACCGGTTAACCGG\n"))
    local_sequences = {str(local): "AACCGGTTAACCGGTT", str(second): "TTGGCCAATTGGCCAA",
                       str(compressed): "TTAACCGGTTAACCGG"}
    attached = []
    for options in [
        {"genome_fasta": str(local)},
        {"genome_fasta": str(second)},
        {"genome_fasta": str(compressed)},
        {"genome_fasta": True},
        {"genome_fasta": True, "genome_fasta_mask": "soft"},
        {"genome_fasta": True, "genome_fasta_mask": "hard"},
    ]:
        with construct(selection, selection.variants[0], sources.root, **options) as genome:
            dna_expected = sources.register(genome, selection.marker)
            if isinstance(options["genome_fasta"], str):
                dna_expected["dna"] = local_sequences[options["genome_fasta"]]
            genome.download()
            genome.index()
            paths = query(genome, dna_expected)
            assert paths[:4] == annotation_paths
            attached.append((genome.to_json(), dna_expected, paths))
    assert len({paths[-1] for _, _, paths in attached}) == 6
    for path, (contents, original) in original_files.items():
        assert path.read_bytes() == contents
        current = path.stat()
        assert (current.st_mode, current.st_ino, current.st_mtime_ns, current.st_ctime_ns) == (
            original.st_mode, original.st_ino, original.st_mtime_ns, original.st_ctime_ns,
        )
    sources.offline = True
    forbid_parsing(monkeypatch)
    modes = {path: path.stat().st_mode for path in sources.root.rglob("*")}
    for path in modes:
        path.chmod(0o555 if path.is_dir() else 0o444)
    before = snapshot(sources.root)
    try:
        for encoded, expected, paths in attached:
            with EnsemblRelease.from_json(encoded) as genome:
                genome.download()
                genome.index()
                assert query(genome, expected) == paths
        assert snapshot(sources.root) == before
    finally:
        for path, mode in modes.items():
            path.chmod(mode)
    for path, (contents, original) in original_files.items():
        assert path.read_bytes() == contents
        current = path.stat()
        # chmod above intentionally changes ctime; attachment operations must
        # preserve original file data, permissions, inode, and modification time.
        assert (current.st_mode, current.st_ino, current.st_mtime_ns) == (
            original.st_mode, original.st_ino, original.st_mtime_ns,
        )


def test_download_only_local_attachments_reopen_without_manifest_changes(sources, tmp_path):
    inputs = [tmp_path / name for name in ("first.fa", "second.fa")]
    for source in inputs:
        source.write_text(">1\nACGT\n")
        with EnsemblRelease(75, genome_fasta=source) as genome:
            genome.download_genome_fasta()
            assert not genome._genome_fasta.index_path.exists()
    sources.offline = True
    before = snapshot(sources.root)
    for source in inputs:
        with EnsemblRelease(75, genome_fasta=source) as genome:
            genome.download_genome_fasta()
    assert snapshot(sources.root) == before


@pytest.mark.parametrize("species, selection_name", [
    ("human", "dated-human-primary"), ("dog", "dated-dog-explicit"),
])
def test_metadata_api_installs_the_same_explicit_dataset(sources, monkeypatch, species, selection_name):
    data = Path(__file__).parent / "data" / "assembly_metadata"
    documents = {
        kind: (data / ("%s-%s.json" % (species, kind))).read_bytes()
        for kind in ("explain", "details", "ftplinks")
    }
    explained = json.loads(documents["explain"])
    genome_id = explained["genome_id"]
    accession = explained["assembly"]["accession_id"]
    endpoints = {
        "/genome/%s/%s" % (accession if kind == "explain" else genome_id, kind): payload
        for kind, payload in documents.items()
    }
    monkeypatch.setattr(
        ensembl_metadata, "urlopen",
        lambda request, **kwargs: BytesIO(
            endpoints[request.full_url.removeprefix(ensembl_metadata.METADATA_API)]
        ),
    )
    selection = next(case for case in SELECTIONS if case.name == selection_name)
    encoded = []
    for identifier in (accession, genome_id):
        with EnsemblRelease.from_ensembl(identifier, species=species, genome_fasta=True) as genome:
            expected = sources.register(genome, selection.marker)
            genome.download()
            genome.index()
            paths = query(genome, expected)
            assert genome.release == selection.release
            encoded.append((genome.to_json(), pickle.dumps(genome)))
    sources.offline = True
    monkeypatch.setattr(
        ensembl_metadata, "urlopen", lambda *a, **kw: pytest.fail("Unexpected metadata lookup")
    )
    forbid_parsing(monkeypatch)
    before = snapshot(sources.root)
    with construct(selection, {"assembly_accession": accession}, sources.root) as genome:
        genome.download()
        genome.index()
        assert query(genome, expected) == paths
    for saved_json, saved_pickle in encoded:
        for genome in (EnsemblRelease.from_json(saved_json), pickle.loads(saved_pickle)):
            with genome:
                genome.download()
                genome.index()
                assert query(genome, expected) == paths
                assert genome.genome_id == genome_id
    assert snapshot(sources.root) == before


def test_conflicting_arguments_do_not_change_an_existing_installation(sources):
    selection = SELECTIONS[1]
    with construct(selection, selection.variants[0], sources.root) as first:
        expected = sources.register(first, selection.marker, selection.alt)
        first.download()
        first.index()
        query(first, expected)
    sources.offline = True
    before = snapshot(sources.root)
    for release, options in [
        (75, {"assembly": "GRCh38"}),
        (75, {"assembly": "GRCh37", "assembly_accession": HUMAN}),
        ("2013_09", {"assembly": "GRCh37.p13", "species": "dog"}),
        ("2019_08", {"species": "dog"}),
        ("2019_08", {"species": "dog_great_dane", "assembly_accession": BASENJI}),
        ("2026_04", {"assembly_accession": HUMAN, "provider": "../outside"}),
    ]:
        with pytest.raises(ValueError):
            EnsemblRelease(release, **options)
    assert snapshot(sources.root) == before
