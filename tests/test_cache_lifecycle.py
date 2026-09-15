import gc
from os.path import exists
from tempfile import TemporaryDirectory
import weakref

from pyensembl import Genome

from .data import (
    MOUSE_ENSMUSG00000017167_PATH,
    MOUSE_ENSMUSG00000017167_PROTEIN_FASTA_PATH,
    MOUSE_ENSMUSG00000017167_TRANSCRIPT_FASTA_PATH,
)


def make_genome(cache_directory_path):
    return Genome(
        reference_name="GRCm38",
        annotation_name="_test_cache_lifecycle",
        gtf_path_or_url=MOUSE_ENSMUSG00000017167_PATH,
        transcript_fasta_paths_or_urls=[
            MOUSE_ENSMUSG00000017167_TRANSCRIPT_FASTA_PATH
        ],
        protein_fasta_paths_or_urls=[
            MOUSE_ENSMUSG00000017167_PROTEIN_FASTA_PATH
        ],
        cache_directory_path=cache_directory_path,
    )


def index_paths(genome):
    return [
        genome.db.local_db_path,
        *genome.transcript_sequences.fasta_dictionary_pickle_paths,
        *genome.protein_sequences.fasta_dictionary_pickle_paths,
    ]


def test_clear_cache_only_clears_in_memory_values():
    with TemporaryDirectory() as cache_directory_path:
        genome = make_genome(cache_directory_path)
        genome.index()
        paths = index_paths(genome)

        genome.gene_by_id("ENSMUSG00000017167")
        genome.transcript_by_id("ENSMUST00000103109")
        genome.exon_by_id("ENSMUSE00000243064")
        first_query_result = genome.db.query(
            select_column_names=["gene_id"],
            filter_column="gene_name",
            filter_value="Cntnap1",
            feature="gene",
        )
        genome.transcript_sequences.get("ENSMUST00000103109")
        genome.protein_sequences.get("ENSMUSP00000099398")

        assert genome._genes
        assert genome._transcripts
        assert genome._exons
        assert all(exists(path) for path in paths)

        genome.clear_cache()

        assert not genome._genes
        assert not genome._transcripts
        assert not genome._exons
        assert genome.transcript_sequences._fasta_dictionary is None
        assert genome.protein_sequences._fasta_dictionary is None
        assert all(exists(path) for path in paths)
        second_query_result = genome.db.query(
            select_column_names=["gene_id"],
            filter_column="gene_name",
            filter_value="Cntnap1",
            feature="gene",
        )
        assert second_query_result == first_query_result
        assert second_query_result is not first_query_result
        genome.db.close()


def test_delete_index_files_preserves_sources_and_allows_reindexing():
    with TemporaryDirectory() as cache_directory_path:
        genome = make_genome(cache_directory_path)
        genome.index()
        paths = index_paths(genome)
        source_paths = [
            MOUSE_ENSMUSG00000017167_PATH,
            MOUSE_ENSMUSG00000017167_TRANSCRIPT_FASTA_PATH,
            MOUSE_ENSMUSG00000017167_PROTEIN_FASTA_PATH,
        ]
        assert all(exists(path) for path in paths)

        genome.delete_index_files()

        assert not any(exists(path) for path in paths)
        assert all(exists(path) for path in source_paths)
        assert genome.db._connection is None

        genome.index()
        assert all(exists(path) for path in paths)
        genome.db.close()


def test_memoized_methods_do_not_keep_genome_alive():
    with TemporaryDirectory() as cache_directory_path:
        genome = make_genome(cache_directory_path)
        genome.index()
        transcript = genome.transcript_by_id("ENSMUST00000103109")
        assert transcript.contains_start_codon
        genome_reference = weakref.ref(genome)
        genome.db.close()

        del transcript
        del genome
        gc.collect()

        assert genome_reference() is None
