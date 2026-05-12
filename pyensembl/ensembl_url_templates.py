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

from .species import Species, find_species_by_name
from .ensembl_versions import check_release_number

ENSEMBL_FTP_SERVER = "https://ftp.ensembl.org"
ENSEMBL_GENOMES_FTP_SERVER = "https://ftp.ensemblgenomes.ebi.ac.uk"
# Back-compat alias; new code should use ENSEMBL_GENOMES_FTP_SERVER.
ENSEMBL_PLANTS_FTP_SERVER = ENSEMBL_GENOMES_FTP_SERVER

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
    return ensembl_release, species.latin_name, reference_name


# GTF annotation file example: Homo_sapiens.GRCh38.gtf.gz
GTF_FILENAME_TEMPLATE = "%(Species)s.%(reference)s.%(release)d.gtf.gz"


def make_gtf_filename(ensembl_release, species):
    """
    Return GTF filename expected on the Ensembl FTP server for a specific
    species/release combination.
    """
    ensembl_release, species_name, reference_name = normalize_release_properties(
        ensembl_release, species
    )
    return GTF_FILENAME_TEMPLATE % {
        "Species": species_name.capitalize(),
        "reference": reference_name,
        "release": ensembl_release,
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
            "species": species.latin_name,
        }
    else:
        if server is None:
            server = ENSEMBL_FTP_SERVER
        subdir = GTF_SUBDIR_TEMPLATE % {
            "release": check_release_number(ensembl_release),
            "species": species.latin_name,
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


def make_fasta_filename(ensembl_release, species, sequence_type, is_plant=None):
    """
    ``is_plant`` is accepted for backward compatibility but ignored; the
    layout is now derived from ``species.ensembl_genomes``.
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
    is_plant=None,
    server=None,
):
    """Construct URL to FASTA file with cDNA transcript or protein sequences.

    ``is_plant`` is accepted for backward compatibility but ignored; routing
    is derived from ``species.ensembl_genomes`` and ``species.division``.
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
