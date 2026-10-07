# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Templates for URLs and paths to specific release, species, and file type
on the Ensembl ftp server (main + Ensembl Genomes).

For example, the human chromosomal DNA sequences for release 78 are in:

    https://ftp.ensembl.org/pub/release-78/fasta/homo_sapiens/dna/
"""

from collections import namedtuple

from .species import Species, find_species_by_name
from .ensembl_versions import check_release_number

ENSEMBL_FTP_SERVER = "https://ftp.ensembl.org"
ENSEMBL_GENOMES_FTP_SERVER = "https://ftp.ensemblgenomes.ebi.ac.uk"
# The new Ensembl platform: dated releases by assembly accession and provider.
ENSEMBL_PLATFORM_FTP_SERVER = "https://ftp.ebi.ac.uk/pub/ensemblorganisms"

# Path layouts:
#   main Ensembl:    /pub/release-N/{gtf,fasta}/{species}/...
#   Ensembl Genomes: /pub/release-N/{division}/{gtf,fasta}/{species}/...
FASTA_SUBDIR_TEMPLATE = "/pub/release-%(release)d/fasta/%(species)s/%(type)s/"
GTF_SUBDIR_TEMPLATE = "/pub/release-%(release)d/gtf/%(species)s/"
GENOMES_FASTA_SUBDIR_TEMPLATE = (
    "/pub/release-%(release)d/%(division)s/fasta/%(species)s/%(type)s/"
)
GENOMES_GTF_SUBDIR_TEMPLATE = (
    "/pub/release-%(release)d/%(division)s/gtf/%(species)s/"
)


def _resolve_species(species):
    if isinstance(species, Species):
        return species
    return find_species_by_name(species)


def normalize_release_properties(ensembl_release, species):
    """
    Make sure a given release is valid, normalize it to be an integer,
    normalize the species name, and get its associated reference.
    """
    ensembl_release = check_release_number(ensembl_release)
    species = _resolve_species(species)
    reference_name = species.which_reference(ensembl_release)
    return ensembl_release, species.ensembl_name(ensembl_release), reference_name


# GTF annotation file example: Homo_sapiens.GRCh38.gtf.gz
GTF_FILENAME_TEMPLATE = "%(Species)s.%(reference)s.%(release)d.gtf.gz"
GTF_FILENAME_TEMPLATE_WITH_PATCHES = (
    "%(Species)s.%(reference)s.%(release)d.chr_patch_hapl_scaff.gtf.gz"
)


def make_gtf_filename(ensembl_release, species):
    """
    Return GTF filename expected on the Ensembl FTP server for a specific
    species/release combination.
    """
    species = _resolve_species(species)
    ensembl_release, species_name, reference_name = normalize_release_properties(
        ensembl_release, species
    )
    # Ensembl split off patch/haplotype annotations in release 82. The full
    # variant exists for GRCh38, GRCm38 (through release 102), and GRCz11
    # (from release 92). Other assemblies, including GRCm39, use the standard
    # file. Resolve the assembly first so transitions retain valid filenames.
    template = GTF_FILENAME_TEMPLATE
    if ensembl_release >= 82 and reference_name in ("GRCh38", "GRCm38", "GRCz11"):
        template = GTF_FILENAME_TEMPLATE_WITH_PATCHES
    # Ensembl 116 republished the Ensembl Genomes 63 GTFs of its
    # non-vertebrates (worm, fly, yeast) under their Ensembl Genomes names.
    filename_release = ensembl_release
    if ensembl_release == 116 and not species.is_vertebrate:
        filename_release = 63
    return template % {
        "Species": species_name.capitalize(),
        "reference": reference_name,
        "release": filename_release,
    }


def make_gtf_url(ensembl_release, species, server=None):
    """
    Returns a fully-qualified URL to the GTF file for ``species`` at
    ``ensembl_release``. Routes through the Ensembl Genomes server when
    ``species.ensembl_genomes`` is True, else through main Ensembl.
    """
    species = _resolve_species(species)
    if species.ensembl_genomes:
        if server is None:
            server = ENSEMBL_GENOMES_FTP_SERVER
        subdir = GENOMES_GTF_SUBDIR_TEMPLATE % {
            "release": check_release_number(ensembl_release),
            "division": species.division,
            "species": species.ensembl_name(check_release_number(ensembl_release)),
        }
    else:
        if server is None:
            server = ENSEMBL_FTP_SERVER
        subdir = GTF_SUBDIR_TEMPLATE % {
            "release": check_release_number(ensembl_release),
            "species": species.ensembl_name(check_release_number(ensembl_release)),
        }
    filename = make_gtf_filename(
        ensembl_release=ensembl_release, species=species
    )
    return server + subdir + filename


# cDNA & protein FASTA file for releases before (and including) Ensembl 75
# example: Homo_sapiens.NCBI36.54.cdna.all.fa.gz
OLD_FASTA_FILENAME_TEMPLATE = (
    "%(Species)s.%(reference)s.%(release)d.%(sequence_type)s.all.fa.gz"
)
OLD_FASTA_FILENAME_TEMPLATE_NCRNA = (
    "%(Species)s.%(reference)s.%(release)d.ncrna.fa.gz"
)
# cDNA & protein FASTA for releases after Ensembl 75 (and all Ensembl Genomes
# releases, which use the modern layout regardless of release number).
NEW_FASTA_FILENAME_TEMPLATE = (
    "%(Species)s.%(reference)s.%(sequence_type)s.all.fa.gz"
)
NEW_FASTA_FILENAME_TEMPLATE_NCRNA = "%(Species)s.%(reference)s.ncrna.fa.gz"


def make_fasta_filename(ensembl_release, species, sequence_type):
    """Filename of a cDNA, ncRNA, or protein FASTA; the layout depends on
    the release and whether the species is on Ensembl Genomes.
    """
    species = _resolve_species(species)
    ensembl_release, species_name, reference_name = normalize_release_properties(
        ensembl_release, species
    )
    use_new_layout = ensembl_release > 75 or species.ensembl_genomes
    if use_new_layout:
        if sequence_type == "ncrna":
            return NEW_FASTA_FILENAME_TEMPLATE_NCRNA % {
                "Species": species_name.capitalize(),
                "reference": reference_name,
            }
        return NEW_FASTA_FILENAME_TEMPLATE % {
            "Species": species_name.capitalize(),
            "reference": reference_name,
            "sequence_type": sequence_type,
        }
    if sequence_type == "ncrna":
        return OLD_FASTA_FILENAME_TEMPLATE_NCRNA % {
            "Species": species_name.capitalize(),
            "reference": reference_name,
            "release": ensembl_release,
        }
    return OLD_FASTA_FILENAME_TEMPLATE % {
        "Species": species_name.capitalize(),
        "reference": reference_name,
        "release": ensembl_release,
        "sequence_type": sequence_type,
    }


def make_fasta_url(
    ensembl_release,
    species,
    sequence_type,
    server=None,
):
    """Construct URL to FASTA file with cDNA transcript or protein sequences.

    Routing is derived from ``species.ensembl_genomes`` and ``species.division``.
    """
    species = _resolve_species(species)
    ensembl_release, species_name, _ = normalize_release_properties(
        ensembl_release, species
    )
    if species.ensembl_genomes:
        if server is None:
            server = ENSEMBL_GENOMES_FTP_SERVER
        subdir = GENOMES_FASTA_SUBDIR_TEMPLATE % {
            "release": ensembl_release,
            "division": species.division,
            "species": species_name,
            "type": sequence_type,
        }
    else:
        if server is None:
            server = ENSEMBL_FTP_SERVER
        subdir = FASTA_SUBDIR_TEMPLATE % {
            "release": ensembl_release,
            "species": species_name,
            "type": sequence_type,
        }
    filename = make_fasta_filename(
        ensembl_release=ensembl_release,
        species=species,
        sequence_type=sequence_type,
    )
    return server + subdir + filename


def make_genome_fasta_url(
    ensembl_release, species, fasta_type="toplevel", mask="none", server=None
):
    """URL for combined reference DNA, including historical/division layouts.

    Toplevel includes patches and haplotypes. Primary assembly excludes them
    and is not available for every species/release. Masked files live in the
    same ``dna`` directory as unmasked files.
    """
    if fasta_type not in ("toplevel", "primary_assembly"):
        raise ValueError("genome_fasta_type must be 'toplevel' or 'primary_assembly'")
    if mask not in ("none", "soft", "hard"):
        raise ValueError("genome_fasta_mask must be 'none', 'soft', or 'hard'")
    species = _resolve_species(species)
    release, species_name, reference = normalize_release_properties(ensembl_release, species)
    prefix = "%s.%s" % (species_name.capitalize(), reference)
    if release <= 75 and not species.ensembl_genomes:
        prefix += ".%d" % release
    sequence_type = {"none": "dna", "soft": "dna_sm", "hard": "dna_rm"}[mask]
    filename = "%s.%s.%s.fa.gz" % (prefix, sequence_type, fasta_type)
    # Reuse the established division routing without treating DNA as cDNA.
    directory = make_fasta_url(release, species, "dna", server=server).rsplit("/", 1)[0]
    return directory + "/" + filename


DatedReleaseUrls = namedtuple(
    "DatedReleaseUrls", ["directory", "gtf", "cdna", "pep", "genome_fasta"]
)
DATED_GENOME_FASTA_MASKS = {"none": "unmasked", "soft": "softmasked", "hard": "hardmasked"}


def make_dated_releases_directory(
    assembly_accession, provider, server=ENSEMBL_PLATFORM_FTP_SERVER
):
    """Directory listing one assembly/provider's annotation dates, e.g.
    .../GCA/000/001/405/29/ensembl for GCA_000001405.29 and "ensembl"."""
    prefix, digits = assembly_accession.split("_")
    number, version = digits.split(".")
    accession_path = "/".join([prefix, number[:3], number[3:6], number[6:9], version])
    return "/".join([server.rstrip("/"), accession_path, provider])


def make_dated_release_urls(
    assembly_accession, provider, annotation_date, include_alt=False,
    genome_fasta_mask="none", server=ENSEMBL_PLATFORM_FTP_SERVER,
):
    """URLs of one assembly/provider/date dataset on the new Ensembl platform.

    For example, GCA_000001405.29 with provider "ensembl" and date "2026_04"
    lives in .../GCA/000/001/405/29/ensembl/2026_04/. cDNA covers every
    transcript biotype, so there is no separate ncRNA FASTA.
    """
    if genome_fasta_mask not in DATED_GENOME_FASTA_MASKS:
        raise ValueError("genome_fasta_mask must be 'none', 'soft', or 'hard'")
    directory = "%s/%s" % (
        make_dated_releases_directory(assembly_accession, provider, server),
        annotation_date,
    )
    geneset = directory + "/geneset/"
    return DatedReleaseUrls(
        directory=directory,
        gtf=geneset + ("genes-including_alt.gtf.gz" if include_alt else "genes.gtf.gz"),
        cdna=geneset + "cdna.fa.bgz",
        pep=geneset + "pep.fa.bgz",
        genome_fasta="%s/genome/%s.fa.bgz" % (
            directory, DATED_GENOME_FASTA_MASKS[genome_fasta_mask]
        ),
    )


def dated_genome_fasta_mask(url):
    """genome_fasta_mask of a dated release's reference DNA URL, or None."""
    for mask, name in DATED_GENOME_FASTA_MASKS.items():
        if url.endswith("/genome/%s.fa.bgz" % name):
            return mask
    return None
