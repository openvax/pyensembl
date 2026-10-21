"""Transcript/CDS source declarations win, including on untouched old indexes."""

import gzip
import json
from pathlib import Path
import sqlite3
import tarfile

import datacache
import pytest

from pyensembl import Genome, EnsemblRelease, Transcript
from pyensembl.annotation_genetic_code import AnnotationGeneticCodes
from pyensembl.database import Database
from .test_legacy_cache_compatibility import snapshot


def write_gtf(path, transcript='', cds='', contig='MT', identifier='ENST1', extra=''):
    attributes = 'gene_id "g"; transcript_id "%s"; transcript_version "2";' % identifier
    text = (
        f'{contig}\ttest\tgene\t1\t9\t.\t+\t.\tgene_id "g"; gene_name "G";\n'
        f'{contig}\ttest\ttranscript\t1\t9\t.\t+\t.\t{attributes} {transcript}\n'
        f'{contig}\ttest\texon\t1\t9\t.\t+\t.\t{attributes} exon_id "e"; exon_number "1";\n'
        f'{contig}\ttest\tCDS\t1\t9\t.\t+\t0\t{attributes} {cds}\n' + extra
    )
    path.write_bytes(gzip.compress(text.encode(), mtime=0) if path.suffix == '.gz' else text.encode())


def open_genome(tmp_path, **kwargs):
    return Genome("GRCh38", "source", gtf_path_or_url=str(tmp_path / "source.gtf"),
                  cache_directory_path=str(tmp_path / "cache"), **kwargs)


@pytest.mark.parametrize("transcript,cds,expected,source", [
    ('', '', 2, 'taxonomy'), ('transl_table "4";', '', 4, 'annotation'),
    ('', 'transl_table "5";', 5, 'annotation'),
    ('translation_table "3";', 'transl_table "3";', 3, 'annotation'),
    ('', 'transl_table "27";', 27, 'annotation'),
    ('', "transl_table '5';", 5, 'annotation'),
    ('', 'transl_table 5;', 5, 'annotation'),
    ('', 'transl_table "5"; transl_table "5";', 5, 'annotation'),
    ('', 'transl_table "";', 2, 'taxonomy'),
    ('', 'note "irrelevant; transl_table 9";', 2, 'taxonomy'),
])
@pytest.mark.parametrize("legacy", [False, True])
def test_source_override_index_or_legacy_stream(tmp_path, monkeypatch, transcript, cds, expected, source, legacy):
    write_gtf(tmp_path / "source.gtf", transcript, cds)
    with open_genome(tmp_path) as genome:
        genome.index()
        if legacy:
            with sqlite3.connect(genome.db.local_db_path) as connection:
                for table in ('gene', 'transcript', 'exon', 'CDS'):
                    for column in ('transl_table', 'translation_table'):
                        connection.execute('ALTER TABLE "%s" DROP COLUMN "%s"' % (table, column))
        else:
            monkeypatch.setattr(AnnotationGeneticCodes, '_read_gtf', lambda *a: pytest.fail('new index should retain declarations'))
        before = snapshot(tmp_path)
        tx = genome.transcript_by_id("ENST1.2")
        assert tx.translation_table == expected
        assert tx.source_translation_table == (expected if source == 'annotation' else None)
        assert tx.genetic_code.source == source
        if source == 'annotation':
            assert tx.genetic_code.source_url == str(tmp_path / 'source.gtf')
        monkeypatch.setattr(AnnotationGeneticCodes, '_read_index', lambda *a: pytest.fail('reuse in-memory result'))
        monkeypatch.setattr(AnnotationGeneticCodes, '_read_gtf', lambda *a: pytest.fail('reuse in-memory result'))
        assert tx.translation_table == expected
        assert snapshot(tmp_path) == before


@pytest.mark.parametrize("transcript,cds,match", [
    ('transl_table "2";', 'transl_table "5";', 'Conflicting'),
    ('', 'transl_table "2"; translation_table "4";', 'Conflicting'),
    ('', 'transl_table "2"; transl_table "4";', 'Conflicting'),
    ('', 'transl_table "999";', 'Unknown NCBI'),
    ('', 'transl_table "not-a-number";', 'Invalid translation table'),
])
@pytest.mark.parametrize("indexed", [False, True])
def test_bad_declarations_do_not_fall_back_to_taxonomy(tmp_path, transcript, cds, match, indexed):
    write_gtf(tmp_path / 'source.gtf', transcript, cds)
    with open_genome(tmp_path) as genome:
        if indexed:
            genome.index()
            tx = genome.transcript_by_id('ENST1')
        else:
            tx = Transcript('ENST1', None, 'MT', 1, 9, '+', None, 'g', genome)
        with pytest.raises(ValueError, match=match):
            _ = tx.translation_table


def test_conflicting_cds_rows_are_not_silently_deduplicated(tmp_path):
    extra = 'MT\ttest\tCDS\t1\t3\t.\t+\t0\tgene_id "g"; transcript_id "ENST1"; transl_table "4";\n'
    write_gtf(tmp_path / 'source.gtf', cds='transl_table "5";', extra=extra)
    with open_genome(tmp_path) as genome:
        genome.index()
        with pytest.raises(ValueError, match='Conflicting'):
            _ = genome.transcript_by_id('ENST1').translation_table


@pytest.mark.parametrize('indexed', [False, True])
def test_source_overrides_reference_on_unplaced_contig(tmp_path, indexed):
    write_gtf(tmp_path / 'source.gtf', cds='transl_table "5";', contig='scaffold_42')
    with open_genome(tmp_path, reference_metadata={'translation_tables': {'nuclear': 4}}) as genome:
        if indexed:
            genome.index()
            tx = genome.transcript_by_id('ENST1')
        else:
            tx = Transcript('ENST1', None, 'scaffold_42', 1, 9, '+', None, 'g', genome)
        result = tx.genetic_code
        assert (result.table_id, result.compartment, result.source) == (5, 'nuclear', 'annotation')


@pytest.mark.parametrize('indexed', [False, True])
@pytest.mark.parametrize('taxon,expected', [(9606, 1), (237561, 12)])
def test_unplaced_transcript_uses_species_default_without_cache_changes(tmp_path, indexed, taxon, expected):
    path = tmp_path / 'source.gtf'
    write_gtf(path, contig='unplaced_42')
    with Genome('custom', 'test', gtf_path_or_url=str(path),
                cache_directory_path=str(tmp_path / 'cache'),
                reference_metadata={'taxon_id': taxon}) as genome:
        if indexed:
            genome.index()
            tx = genome.transcript_by_id('ENST1')
        else:
            tx = Transcript('ENST1', None, 'unplaced_42', 1, 9, '+', None, 'g', genome)
        before = snapshot(tmp_path)
        assert tx.source_translation_table is None
        assert tx.translation_table == expected
        assert tx.genetic_code.compartment == 'nuclear'
        assert tx.genetic_code.source == 'taxonomy'
        assert snapshot(tmp_path) == before


@pytest.mark.parametrize('compressed', [False, True])
def test_sparse_legacy_read_and_source_change_invalidation(tmp_path, compressed):
    path = tmp_path / ('source.gtf.gz' if compressed else 'source.gtf')
    write_gtf(path, cds='transl_table "5";', identifier='mRNA-Solyc04g011850.1.1')
    reader = AnnotationGeneticCodes(path, path, str(tmp_path))
    before = snapshot(tmp_path)
    assert reader.get('mRNA-Solyc04g011850.1.1') == 5
    assert reader.get('mRNA-Solyc04g011850.1') is None
    assert snapshot(tmp_path) == before
    write_gtf(path, cds='transl_table "4";', identifier='mRNA-Solyc04g011850.1.1')
    assert reader.get('mRNA-Solyc04g011850.1.1') == 4


def test_changed_remote_source_does_not_reuse_stale_index_metadata(tmp_path):
    path = tmp_path / 'source.gtf'
    write_gtf(path, cds='transl_table "5";')
    options = dict(gtf_path_or_url='https://example.invalid/source.gtf', cache_directory_path=str(tmp_path))
    with Genome('GRCh38', 'remote', **options) as genome:
        genome.index()
        tx = genome.transcript_by_id('ENST1')
        assert tx.translation_table == 5
        write_gtf(path, cds='transl_table "4";')
        assert tx.translation_table == 4
    with Genome('GRCh38', 'remote', **options) as reopened:
        assert reopened.transcript_by_id('ENST1').translation_table == 4


def test_restricted_custom_index_does_not_hide_cds_declaration(tmp_path):
    write_gtf(tmp_path / 'source.gtf', cds='transl_table "5";')
    with open_genome(tmp_path) as genome:
        genome.db.restrict_gtf_features = {'gene', 'transcript'}
        genome.index()
        assert genome.transcript_by_id('ENST1').translation_table == 5


def test_missing_annotation_does_not_download_or_assume_no_override(tmp_path, monkeypatch):
    monkeypatch.setattr(datacache, 'fetch_file', lambda *a, **k: pytest.fail('network'))
    genome = Genome('GRCh38', 'missing', gtf_path_or_url='https://example.invalid/genes.gtf',
                    cache_directory_path=str(tmp_path))
    tx = Transcript('t', None, 'MT', 1, 9, '+', None, 'g', genome)
    with pytest.raises(FileNotFoundError, match='installed GTF is missing'):
        _ = tx.translation_table
    assert not list(tmp_path.iterdir())
    no_gtf = Genome('GRCh38', 'metadata-only')
    assert Transcript('t', None, 'MT', 1, 9, '+', None, 'g', no_gtf).translation_table == 2


def test_published_230_indexes_stay_read_only_and_byte_identical(tmp_path, monkeypatch):
    data = Path(__file__).parent / 'data/tomato_sl4'
    with tarfile.open(data / 'pyensembl-2.30.0-cache.tar.gz') as archive:
        archive.extractall(tmp_path, filter='data')
    root = tmp_path / 'cache'
    monkeypatch.setenv('PYENSEMBL_CACHE_DIR', str(root))
    monkeypatch.setattr(datacache, 'fetch_file', lambda *a, **k: pytest.fail('network'))
    monkeypatch.setattr(Database, 'create', lambda *a, **k: pytest.fail('reindex'))
    modes = {path: path.stat().st_mode for path in root.rglob('*')}
    for path in modes:
        path.chmod(0o555 if path.is_dir() else 0o444)
    before = snapshot(root)
    try:
        for release in (62, '2018_10'):
            with EnsemblRelease(release, species='tomato') as genome:
                tx = genome.transcripts()[0]
                assert tx.source_translation_table is None
                assert tx.translation_table == 1
                assert not genome.db.column_exists('transcript', 'transl_table')
        assert snapshot(root) == before
    finally:
        for path, mode in modes.items():
            path.chmod(mode)


def test_actual_2301_index_finds_previously_unindexed_source_override(tmp_path, monkeypatch):
    data = Path(__file__).parent / 'data/genetic_codes'
    with tarfile.open(data / 'pyensembl-2.30.1-cache.tar.gz') as archive:
        archive.extractall(tmp_path, filter='data')
    root = tmp_path / 'cache'
    options = json.loads((data / 'pyensembl-2.30.1.json').read_text())
    options['cache_directory_path'] = str(root)
    monkeypatch.setattr(datacache, 'fetch_file', lambda *a, **k: pytest.fail('network'))
    monkeypatch.setattr(Database, 'create', lambda *a, **k: pytest.fail('reindex'))
    modes = {path: path.stat().st_mode for path in root.rglob('*')}
    for path in modes:
        path.chmod(0o555 if path.is_dir() else 0o444)
    before = snapshot(root)
    try:
        with Genome.from_dict(options) as genome:
            assert not genome.db.column_exists('CDS', 'transl_table')
            tx = genome.transcript_by_id('ENST1.2')
            assert genome.genetic_code('MT').table_id == 2
            assert tx.source_translation_table == tx.translation_table == 5
            assert tx.genetic_code.source == 'annotation'
        assert snapshot(root) == before
    finally:
        for path, mode in modes.items():
            path.chmod(mode)
