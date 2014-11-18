from pyensembl.locus import normalize_chromosome, Locus

from nose.tools import assert_raises


def test_normalize_chromosome():
    assert normalize_chromosome("X") == "X"
    assert normalize_chromosome("chrX") == "X"
    assert normalize_chromosome("x") == "X"

    assert normalize_chromosome(1) == "1"
    assert normalize_chromosome("1") == "1"
    assert normalize_chromosome("chr1") == "1"

    assert normalize_chromosome("chrM") == "MT"
    assert normalize_chromosome("chrMT") == "MT"
    assert normalize_chromosome("M") == "MT"
    assert normalize_chromosome("MT") == "MT"

    with assert_raises(TypeError):
        normalize_chromosome({"a":"b"})

    with assert_raises(TypeError):
        normalize_chromosome([])

    with assert_raises(TypeError):
        normalize_chromosome(None)

    with assert_raises(ValueError):
        normalize_chromosome("")

    with assert_raises(ValueError):
        normalize_chromosome(0)

def test_locus_overlaps():
    locus = Locus("1", 10, 20, "+")
    assert locus.overlaps("1", 10, 20, "+")
    assert locus.overlaps("1", 10, 20)
    assert locus.overlaps("1", 5, 30)
    assert locus.overlaps("1", 15, 16)
    assert locus.overlaps("1", 15, 30)
    assert locus.overlaps("1", 5, 15)
    assert locus.overlaps("1", 10, 10)
    assert locus.overlaps("1", 20, 20)
    # before start
    assert not locus.overlaps(1,9)
    # after end
    assert not locus.overlaps(21, 30)
    # wrong contig
    assert not locus.overlaps("2", 10, 20)
    # wrong strand
    assert not locus.overlaps("1", 10, 20, "-")

def test_locus_overlaps():
    locus = Locus("1", 10, 20, "+")
    assert locus.overlaps("1", 10, 20, "+")
    assert locus.overlaps("1", 10, 20)
    assert locus.overlaps("1", 15, 16)
    assert locus.overlaps("1", 10, 10)
    assert locus.overlaps("1", 20, 20)
    assert locus.overlaps("1", 5, 30)
    assert locus.overlaps("1", 15, 30)
    assert locus.overlaps("1", 5, 15)

    # before start
    assert not locus.overlaps("1", 1,9)

    # after end
    assert not locus.overlaps("1", 21, 30)

    # wrong contig
    assert not locus.overlaps("2", 10, 20)

    # wrong strand
    assert not locus.overlaps("1", 10, 20, "-")


def test_locus_contains():
    locus = Locus("1", 10, 20, "+")
    assert locus.contains("1", 10, 20, "+")
    assert locus.contains("1", 10, 20)
    assert locus.contains("1", 15, 16)
    assert locus.contains("1", 10, 10)
    assert locus.contains("1", 20, 20)


    # before start and after end
    assert not locus.contains("1", 5, 30)

    # before start
    assert not locus.contains("1", 1,9)
    assert not locus.contains("1", 5, 15)

    # after end
    assert not locus.contains("1", 21, 30)
    assert not locus.contains("1", 15, 30)

    # wrong contig
    assert not locus.contains("2", 10, 20)

    # wrong strand
    assert not locus.contains("1", 10, 20, "-")

