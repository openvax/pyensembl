from __future__ import absolute_import

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
        normalize_chromosome({"a": "b"})

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
    assert not locus.overlaps(1, 9, 9)
    # after end
    assert not locus.overlaps(21, 30, 30)
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
    assert not locus.contains("1", 1, 9)
    assert not locus.contains("1", 5, 15)

    # after end
    assert not locus.contains("1", 21, 30)
    assert not locus.contains("1", 15, 30)

    # wrong contig
    assert not locus.contains("2", 10, 20)

    # wrong strand
    assert not locus.contains("1", 10, 20, "-")

def test_position_offset():
    forward_locus = Locus("1", 10, 20, "+")
    assert forward_locus.offset(10) == 0
    assert forward_locus.offset(15) == 5
    assert forward_locus.offset(19) == 9
    assert forward_locus.offset(20) == 10

    negative_locus = Locus("1", 10, 20, "-")
    assert negative_locus.offset(10) == 10
    assert negative_locus.offset(15) == 5
    assert negative_locus.offset(19) == 1
    assert negative_locus.offset(20) == 0

    # don't allow negative offsets
    with assert_raises(ValueError):
        forward_locus.offset(9)

    # don't allow negative offsets
    with assert_raises(ValueError):
        negative_locus.offset(9)

    # don't allow offset past the end of the locus
    with assert_raises(ValueError):
        forward_locus.offset(21)

    # don't allow offset past the end of the locus
    with assert_raises(ValueError):
        negative_locus.offset(21)


def test_range_offset():
    forward_locus = Locus("1", 10, 20, "+")
    assert forward_locus.offset_range(10, 20) == (0, 10)
    assert forward_locus.offset_range(11, 14) == (1, 4)
    assert forward_locus.offset_range(20, 20) == (10, 10)

    negative_locus = Locus("1", 10, 20, "-")
    assert negative_locus.offset_range(10, 20) == (0, 10)
    assert negative_locus.offset_range(11, 14) == (6, 9)
    assert negative_locus.offset_range(20, 20) == (0, 0)

    # start shouldn't be larger than end
    with assert_raises(AssertionError):
        forward_locus.offset_range(21, 20)

    # start shouldn't be larger than end
    with assert_raises(AssertionError):
        negative_locus.offset_range(21, 20)

    # don't allow negative offsets
    with assert_raises(ValueError):
        forward_locus.offset_range(9, 10)

    # don't allow negative offsets
    with assert_raises(ValueError):
        forward_locus.offset_range(9, 10)

    # don't allow negative offsets
    with assert_raises(ValueError):
        negative_locus.offset_range(9, 10)

def test_locus_distance():
    locus_chr1_10_20_pos = Locus("1", 10, 20, "+")
    locus_chr1_21_25_pos = Locus("1", 21, 25, "+")
    locus_chr2_21_25_pos = Locus("2", 21, 25, "+")
    locus_chr1_21_25_neg = Locus("1", 21, 25, "-")
    assert locus_chr1_10_20_pos.distance_to_locus(locus_chr1_21_25_pos) == 1
    assert locus_chr1_21_25_pos.distance_to_locus(locus_chr1_10_20_pos) == 1
    inf = float("inf")
    assert locus_chr1_10_20_pos.distance_to_locus(locus_chr2_21_25_pos) == inf
    assert locus_chr1_10_20_pos.distance_to_locus(locus_chr1_21_25_neg) == inf
