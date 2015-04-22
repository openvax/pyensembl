from __future__ import print_function, absolute_import

from pyensembl import ensembl_grch38 as ensembl

from tinytimer import benchmark

contigs = [str(i + 1) for i in range(22)] + ["X", "Y"]

def make_repeat_lookup_fn(lookup_fn, n_positions):
    """
    Make a thunk which calls the lookup_fn at a number of loci
    for each human chromosome (excluding MT).
    """
    def repeat_lookup_fn():
        for contig in contigs:
            for position in [10 ** 6 + i * 10 ** 6 for i in range(n_positions)]:
                lookup_fn(contig, position)
    return repeat_lookup_fn

def run_benchmark(lookup_fn, n_positions_per_contig=20, time_limit=60.0):
    """
    Take a lookup functions (such as EnsemblRelease.genes_at_locus) and
    time how long it takes across multiple loci.
    """
    repeat_lookup_fn = make_repeat_lookup_fn(lookup_fn, n_positions_per_contig)
    n_loci = n_positions_per_contig * len(contigs)
    name = lookup_fn.__name__
    average_time = benchmark(
        repeat_lookup_fn,
        name="%s for %d loci" % (name, n_loci))
    print("-- %s : %0.4fs" % (name, average_time))
    assert average_time < time_limit, \
        "%s took too long for %s loci: %0.4fs" % (name, n_loci, average_time)
    return average_time

def test_timing_genes_at_locus():
    run_benchmark(ensembl.genes_at_locus)

def test_timing_transcripts_at_locus():
    run_benchmark(ensembl.transcripts_at_locus)

def test_timing_exons_at_locus():
    run_benchmark(ensembl.exons_at_locus)

def test_timing_transcript_sequences_at_locus():
    def transcript_sequences_at_locus(contig, position):
        sequences = []
        for transcript in ensembl.transcripts_at_locus(contig, position):
            sequences.append(transcript.sequence)
        return sequences
    run_benchmark(transcript_sequences_at_locus)

def test_timing_transcript_coding_sequences_at_locus():
    def transcript_coding_sequences_at_locus(contig, position):
        sequences = []
        for transcript in ensembl.transcripts_at_locus(contig, position):
            if transcript.sequence and transcript.complete:
                sequences.append(transcript.coding_sequence)
        return sequences
    run_benchmark(transcript_coding_sequences_at_locus)

def run_all_benchmarks():
    import types

    # run all local test functions to see their timings printed
    global_variables = globals()
    for variable_name in global_variables:
        if "test_" in variable_name:
            f = global_variables[variable_name]
            if isinstance(f, types.FunctionType):
                f()

if __name__ == "__main__":
    run_all_benchmarks()
