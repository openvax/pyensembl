"""Assembly release ranges agree with Ensembl's archive (#396).

Expected filenames were read from the archive's GTF and FASTA directories.
Set PYENSEMBL_NETWORK_TESTS=1 to recheck every assembly boundary against the
live FTP servers, e.g. after raising MAX_ENSEMBL_RELEASE.
"""

from concurrent.futures import ThreadPoolExecutor
import os
import urllib.error
import urllib.request

import pytest

from pyensembl import EnsemblRelease
from pyensembl.species import Species


@pytest.mark.parametrize(
    "species,release,gtf,cdna",
    [
        # The reported case: release 97 is GRCg6a, not Gallus_gallus-5.0.
        ("chicken", 94, "Gallus_gallus.Gallus_gallus-5.0.94.gtf.gz",
         "Gallus_gallus.Gallus_gallus-5.0.cdna.all.fa.gz"),
        ("chicken", 97, "Gallus_gallus.GRCg6a.97.gtf.gz",
         "Gallus_gallus.GRCg6a.cdna.all.fa.gz"),
        ("chicken", 107, "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.107.gtf.gz",
         "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa.gz"),
        ("cat", 114, "Felis_catus.F.catus_Fca126_mat1.0.114.gtf.gz",
         "Felis_catus.F.catus_Fca126_mat1.0.cdna.all.fa.gz"),
        ("macaca_mulatta", 97, "Macaca_mulatta.Mmul_8.0.1.97.gtf.gz",
         "Macaca_mulatta.Mmul_8.0.1.cdna.all.fa.gz"),
        ("pig", 89, "Sus_scrofa.Sscrofa10.2.89.gtf.gz",
         "Sus_scrofa.Sscrofa10.2.cdna.all.fa.gz"),
        ("drosophila_melanogaster", 111, "Drosophila_melanogaster.BDGP6.46.111.gtf.gz",
         "Drosophila_melanogaster.BDGP6.46.cdna.all.fa.gz"),
        ("naked_mole_rat", 110,
         "Heterocephalus_glaber_female.Naked_mole-rat_maternal.110.gtf.gz",
         "Heterocephalus_glaber_female.Naked_mole-rat_maternal.cdna.all.fa.gz"),
    ],
)
def test_assembly_transitions_name_archived_files(species, release, gtf, cdna):
    genome = EnsemblRelease(release, species=species)
    assert os.path.basename(genome.gtf_url) == gtf
    assert os.path.basename(genome.transcript_fasta_urls[0]) == cdna


@pytest.mark.parametrize(
    "species,release",
    [
        ("dog", 100),  # Later dog releases use canis_lupus_familiaris.
        ("meriones_unguiculatus", 95),  # Added in release 96.
        ("mus_musculus_balbcj", 87),  # 87-91 kept release 86's filenames.
    ],
)
def test_releases_missing_from_the_archive_are_rejected(species, release):
    with pytest.raises(ValueError, match="No genome for"):
        EnsemblRelease(release, species=species)


def _boundary_urls():
    for species in Species._latin_names_to_species.values():
        for start, end in species.reference_assemblies.values():
            for release in sorted({start, end}):
                genome = EnsemblRelease(release, species=species)
                yield genome.gtf_url
                yield genome.transcript_fasta_urls[0]
                yield genome.protein_fasta_urls[0]


def _missing(url):
    for _ in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60):
                return None
        except urllib.error.HTTPError as error:
            return "%s: HTTP %d" % (url, error.code)
        except OSError as error:
            problem = "%s: %s" % (url, error)
    return problem


@pytest.mark.skipif(
    not os.environ.get("PYENSEMBL_NETWORK_TESTS"),
    reason="set PYENSEMBL_NETWORK_TESTS=1 to check the live Ensembl archive",
)
def test_every_assembly_boundary_exists_in_the_archive():
    with ThreadPoolExecutor(6) as pool:
        problems = [problem for problem in pool.map(_missing, _boundary_urls()) if problem]
    assert problems == []
