# Copyright (c) 2015. Mount Sinai School of Medicine
#
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
Templates for URLs and paths to specific relase, species, and file type
on the Ensembl ftp server.

For example, the human chromosomal DNA sequences for release 78 are in:

    ftp://ftp.ensembl.org/pub/release-78/fasta/homo_sapiens/dna/

"""
from __future__ import print_function, division, absolute_import

from os.path import join

from six.moves import urllib_parse

from .species import Species, find_species_by_name
from .ensembl_release_versions import check_release_number

ENSEMBL_FTP_SERVER = "ftp://ftp.ensembl.org"

# Example directories
# FASTA files: /pub/release-78/fasta/homo_sapiens/
# GTF annotation files: /pub/release-78/gtf/homo_sapiens/
SPECIES_SUBDIR_TEMPLATE = "/pub/release-%(release)d/%(filetype)s/%(species)s/"


def _species_subdir(
        ensembl_release,
        species="homo_sapiens",
        filetype="gtf",
        server=ENSEMBL_FTP_SERVER):
    """
    Assume ensembl_release has already been normalize by calling function
    but species might be either a common name or latin name.
    """
    return SPECIES_SUBDIR_TEMPLATE % {
        "release": ensembl_release,
        "filetype": filetype,
        "species": species,
    }


def normalize_release_properties(ensembl_release, species):
    """
    Make sure a given release is valid, normalize it to be an integer,
    normalize the species name, and get its associated reference.
    """
    ensembl_release = check_release_number(ensembl_release)
    if not isinstance(species, Species):
        species = find_species_by_name(species)
    reference_name = species.which_reference(ensembl_release)
    return ensembl_release, species.latin_name, reference_name

# GTF annotation file example: Homo_sapiens.GTCh38.gtf.gz
GTF_FILENAME_TEMPLATE = "%(Species)s.%(reference)s.%(release)d.gtf.gz"


def make_gtf_filename(ensembl_release, species):
    """
    Return GTF filename expect on Ensembl FTP server for a specific
    species/release combination
    """
    ensembl_release, species, reference_name = normalize_release_properties(
        ensembl_release, species)
    return GTF_FILENAME_TEMPLATE % {
        "Species": species.capitalize(),
        "reference": reference_name,
        "release": ensembl_release,
    }

def make_gtf_url(ensembl_release, species, server=ENSEMBL_FTP_SERVER):
    """
    Returns a URL and a filename, which can be joined together.
    """
    ensembl_release, species, _ = \
        normalize_release_properties(ensembl_release, species)
    subdir = _species_subdir(
        ensembl_release,
        species=species,
        filetype="gtf",
        server=server)
    url_subdir = urllib_parse.urljoin(server, subdir)
    filename = make_gtf_filename(
        ensembl_release=ensembl_release,
        species=species)
    return join(url_subdir, filename)

# DNA fasta file example: Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
FASTA_DNA_CHROMOSOME_FILENAME_TEMPLATE = \
    "%(Species)s.%(reference)s.%(release)d.%(sequence_type)s.chromosome.%(contig)s.fa.gz"


def make_fasta_dna_filename(ensembl_release, species, contig):
    ensembl_release, species, reference_name = \
        normalize_release_properties(ensembl_release, species)
    return FASTA_DNA_CHROMOSOME_FILENAME_TEMPLATE % {
        "Species": species.capitalize(),
        "reference": reference_name,
        "release": ensembl_release,
        "sequence_type": "dna",
        "contig": contig
    }

def make_fasta_dna_url(
        ensembl_release,
        species,
        contig,
        server=ENSEMBL_FTP_SERVER):
    """
    Construct URL to FASTA file with full sequence of a particular chromosome.
    Returns server_url/subdir and filename as tuple result.
    """
    ensembl_release, species, _ = \
        normalize_release_properties(ensembl_release, species)
    subdir = _species_subdir(
        ensembl_release,
        species=species,
        filetype="fasta",
        server=server,)
    server_subdir = urllib_parse.urljoin(server, subdir)
    server_sequence_subdir = join(server_subdir, "dna")
    filename = make_fasta_dna_filename(
        ensembl_release=ensembl_release,
        species=species,
        contig=contig)
    return join(server_sequence_subdir, filename)


# cDNA & protein FASTA file for releases before (and including) Ensembl 75
# example: Homo_sapiens.NCBI36.54.cdna.all.fa.gz
OLD_FASTA_FILENAME_TEMPLATE = \
    "%(Species)s.%(reference)s.%(release)d.%(sequence_type)s.all.fa.gz"

# ncRNA FASTA file for releases before (and including) Ensembl 75
# example: Homo_sapiens.NCBI36.54.ncrna.fa.gz

OLD_FASTA_FILENAME_TEMPLATE_NCRNA = \
    "%(Species)s.%(reference)s.%(release)d.ncrna.fa.gz"

# cDNA & protein FASTA file for releases after Ensembl 75
# example: Homo_sapiens.GRCh37.cdna.all.fa.gz
NEW_FASTA_FILENAME_TEMPLATE = \
    "%(Species)s.%(reference)s.%(sequence_type)s.all.fa.gz"

# ncRNA FASTA file for releases after Ensembl 75
# example: Homo_sapiens.GRCh37.ncrna.fa.gz
NEW_FASTA_FILENAME_TEMPLATE_NCRNA = \
    "%(Species)s.%(reference)s.ncrna.fa.gz"


def make_fasta_filename(ensembl_release, species, sequence_type):
    ensembl_release, species, reference_name = \
        normalize_release_properties(ensembl_release, species)
    if ensembl_release <= 75:
        if sequence_type == 'ncrna':
            return OLD_FASTA_FILENAME_TEMPLATE_NCRNA % {
                "Species": species.capitalize(),
                "reference": reference_name,
                "release": ensembl_release
            }
        else:
            return OLD_FASTA_FILENAME_TEMPLATE % {
                "Species": species.capitalize(),
                "reference": reference_name,
                "release": ensembl_release,
                "sequence_type": sequence_type,
            }
    else:
        if sequence_type == 'ncrna':
            return NEW_FASTA_FILENAME_TEMPLATE_NCRNA % {
                "Species": species.capitalize(),
                "reference": reference_name
            }
        else:
            return NEW_FASTA_FILENAME_TEMPLATE % {
                "Species": species.capitalize(),
                "reference": reference_name,
                "sequence_type": sequence_type,
            }

def make_fasta_url(
        ensembl_release,
        species,
        sequence_type,
        server=ENSEMBL_FTP_SERVER):
    """Construct URL to FASTA file with cDNA transcript or protein sequences

    Parameter examples:
        ensembl_release = 75
        species = "Homo_sapiens"
        sequence_type = "cdna" (other option: "pep")
    """
    ensembl_release, species, reference_name = normalize_release_properties(
        ensembl_release, species)
    subdir = _species_subdir(
        ensembl_release,
        species=species,
        filetype="fasta",
        server=server)
    server_subdir = urllib_parse.urljoin(server, subdir)
    server_sequence_subdir = join(server_subdir, sequence_type)
    filename = make_fasta_filename(
        ensembl_release=ensembl_release,
        species=species,
        sequence_type=sequence_type)
    return join(server_sequence_subdir, filename)
