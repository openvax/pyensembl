"""Unit tests for pyensembl.fasta."""

import gzip
import logging

import pytest

from pyensembl.fasta import (
    FastaParser,
    _parse_header_id,
    first_fasta_id,
    parse_fasta_dictionary,
)


def test_parse_header_id_validates_bytes_and_empty_line():
    with pytest.raises(TypeError, match="Expected header line to be of type"):
        _parse_header_id(">ENSP00000123456.1")
    with pytest.raises(ValueError, match="No identifier on FASTA line"):
        _parse_header_id(b">")
    with pytest.raises(ValueError, match="No identifier on FASTA line"):
        _parse_header_id(b"")
    assert _parse_header_id(b"> leading-space") == ""
    assert _parse_header_id(b">|pipe-first") == ""


def test_fasta_parser_ensembl_quirks(tmp_path, caplog):
    fasta = tmp_path / "quirks.fa"
    fasta.write_text(
        "pre-header noise line\n"
        "; semicolon comment before first record\n"
        ">ENSP000001.2 pep:protein_coding\n"
        "; comment inside sequence\n"
        "MSEQ\n"
        "PART2\n"
        "> empty_id_leading_space\n"
        "IGNORED\n"
        ">|empty_id_pipe\n"
        "IGNORED2\n"
        ">ENSP000002.1\n"
        ">ENSP000003.4|ENST000003.1|GENE\n"
        "MFINAL\n"
    )
    with caplog.at_level(logging.WARNING, logger="pyensembl.fasta"):
        parsed = parse_fasta_dictionary(str(fasta))

    assert parsed == {
        "ENSP000001.2": "MSEQPART2",
        "ENSP000003.4": "MFINAL",
    }
    messages = [record.getMessage() for record in caplog.records]
    assert any("Unable to parse ID from header line" in m for m in messages)
    assert any("No sequence data for 'ENSP000002.1'" in m for m in messages)


def test_fasta_parser_bare_gt_raises_no_identifier(tmp_path):
    fasta = tmp_path / "empty_gt.fa"
    fasta.write_text(">\nACGT\n")
    with pytest.raises(ValueError, match="No identifier on FASTA line"):
        parse_fasta_dictionary(str(fasta))

    fasta_spaces = tmp_path / "spaces_gt.fa"
    fasta_spaces.write_text(">   \nACGT\n")
    with pytest.raises(ValueError, match="No identifier on FASTA line"):
        parse_fasta_dictionary(str(fasta_spaces))


@pytest.mark.parametrize("suffix", [".fa", ".fa.gz", ".fa.gzip"])
def test_fasta_parser_and_first_fasta_id_plain_and_gzipped(tmp_path, suffix):
    content = b">ENSP00000493376.2|ENST00000641515.2|GENE\nMSEQ\n>ENSP00000999999.1 desc\nMNEXT\n"
    path = tmp_path / ("proteins" + suffix)
    if suffix.endswith(("gz", "gzip")):
        path.write_bytes(gzip.compress(content))
    else:
        path.write_bytes(content)

    assert first_fasta_id(str(path)) == "ENSP00000493376.2"
    assert FastaParser().read_file(str(path)) == {
        "ENSP00000493376.2": "MSEQ",
        "ENSP00000999999.1": "MNEXT",
    }


def test_first_fasta_id_empty_file_returns_none(tmp_path):
    empty = tmp_path / "empty.fa"
    empty.write_text("")
    assert first_fasta_id(str(empty)) is None
