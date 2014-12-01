from pyensembl import EnsemblRelease

releases = [
    # last release for GRCh36/hg18
    EnsemblRelease(54),
    # last release for GRCh37/hg19
    EnsemblRelease(75),
    # most recent release for GRCh38
    EnsemblRelease(77)
]

contigs = range(1,23) + ["X", "Y", "M"]
