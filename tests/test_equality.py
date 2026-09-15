"""Equality and dictionary-key behavior for public classes and subclasses."""

import pytest

from pyensembl import EnsemblRelease, Gene, Genome, Locus, SequenceData, Transcript
from pyensembl.database import Database
from pyensembl.download_cache import DownloadCache
from pyensembl.protein import Protein
from pyensembl.species import Species


@pytest.fixture(
    params=[
        Gene, Transcript, Protein, Genome, EnsemblRelease, Species,
        SequenceData, Database, DownloadCache,
    ],
    ids=lambda cls: cls.__name__,
)
def equality_factory(request, tmp_path, monkeypatch):
    base_type = request.param
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    gtf = tmp_path / "annotations.gtf"
    gtf.touch()
    for name in ("same", "different"):
        (tmp_path / (name + ".fa")).write_text(">sequence\nACGT\n")
    genome = Genome(
        reference_name="test",
        annotation_name="test",
        gtf_path_or_url=str(gtf),
        cache_directory_path=str(tmp_path),
    )
    locus_fields = dict(
        contig="1", start=10, end=20, strand="+",
        biotype="protein_coding", genome=genome,
    )

    def make(cls, different=False):
        name = "different" if different else "same"
        arguments = {
            Gene: dict(gene_id=name, gene_name=name, **locus_fields),
            Transcript: dict(
                transcript_id=name, transcript_name=name, gene_id="gene",
                **locus_fields,
            ),
            Protein: dict(protein_id=name),
            Genome: dict(
                reference_name="test", annotation_name=name,
                cache_directory_path=str(tmp_path),
            ),
            EnsemblRelease: dict(release=82 if different else 81),
            Species: dict(latin_name=name),
            SequenceData: dict(fasta_paths=[str(tmp_path / (name + ".fa"))]),
            Database: dict(gtf_path=str(tmp_path / (name + ".gtf"))),
            DownloadCache: dict(
                reference_name="test", annotation_name=name,
                cache_directory_path=str(tmp_path),
            ),
        }
        return cls(**arguments[base_type])

    yield base_type, make
    genome.close()


@pytest.mark.parametrize("use_subclass", [False, True], ids=["base", "subclass"])
def test_equality_and_hashing(equality_factory, use_subclass):
    base_type, make = equality_factory
    cls = type("Custom" + base_type.__name__, (base_type,), {}) if use_subclass else base_type
    first = make(cls)
    equal = make(cls)
    different = make(cls, different=True)

    assert first == first
    assert first == equal
    assert equal == first
    assert not first != equal
    assert first != different
    assert different != first
    assert hash(first) == hash(equal)
    assert {first: "cached"}[equal] == "cached"
    assert len({first, equal, different}) == 2


def test_equality_distinguishes_concrete_types(equality_factory):
    base_type, make = equality_factory
    child = type("Custom" + base_type.__name__, (base_type,), {})
    sibling = type("Other" + base_type.__name__, (base_type,), {})
    objects = [make(base_type), make(child), make(sibling)]

    for index, first in enumerate(objects):
        for second in objects[index + 1:] + [None, object()]:
            assert first != second
            assert second != first
            assert not first == second
            assert not second == first
        if isinstance(first, Locus):
            # Returning NotImplemented here would fall through to Locus's
            # coordinate equality and accidentally equate distinct types.
            locus = Locus("1", 10, 20, "+")
            assert first != locus
            assert locus != first

    assert len(set(objects)) == 3
