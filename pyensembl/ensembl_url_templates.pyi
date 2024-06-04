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

    https://ftp.ensembl.org/pub/release-78/fasta/homo_sapiens/dna/

"""

from typing import Literal, Tuple, Union, TYPE_CHECKING
if TYPE_CHECKING:
    from .species import Species

ENSEMBL_FTP_SERVER: str = "https://ftp.ensembl.org"
ENSEMBL_PLANTS_FTP_SERVER: str = "https://ftp.ensemblgenomes.ebi.ac.uk/"

FASTA_SUBDIR_TEMPLATE: str = "/pub/release-%(release)d/fasta/%(species)s/%(type)s/"
PLANTS_FASTA_SUBDIR_TEMPLATE: str = (
    "/pub/release-%(release)d/plants/fasta/%(species)s/%(type)s/"
)
GTF_SUBDIR_TEMPLATE: str = "/pub/release-%(release)d/gtf/%(species)s/"
PLANTS_GTF_SUBDIR_TEMPLATE: str = "/pub/release-%(release)d/plants/gtf/%(species)s/"

lPlants: Tuple[str] = ("arabidopsis_thaliana", "arabidopsis")

def normalize_release_properties(
    ensembl_release: Union[str, int], species: Union[str, Species]
) -> Tuple[int, str, str]: ...

GTF_FILENAME_TEMPLATE: str = "%(Species)s.%(reference)s.%(release)d.gtf.gz"

def make_gtf_filename(
    ensembl_release: Union[str, int], species: Union[str, Species]
) -> str: ...
def make_gtf_url(
    ensembl_release: Union[str, int],
    species: Union[str, Species],
    server: str = ENSEMBL_FTP_SERVER,
    gtf_subdir=GTF_SUBDIR_TEMPLATE,
) -> str: ...

OLD_FASTA_FILENAME_TEMPLATE: str = (
    "%(Species)s.%(reference)s.%(release)d.%(sequence_type)s.all.fa.gz"
)

OLD_FASTA_FILENAME_TEMPLATE_NCRNA: str = (
    "%(Species)s.%(reference)s.%(release)d.ncrna.fa.gz"
)

NEW_FASTA_FILENAME_TEMPLATE: str = (
    "%(Species)s.%(reference)s.%(sequence_type)s.all.fa.gz"
)

NEW_FASTA_FILENAME_TEMPLATE_NCRNA: str = "%(Species)s.%(reference)s.ncrna.fa.gz"

def make_fasta_filename(
    ensembl_release: Union[str, int],
    species: Union[str, Species],
    sequence_type: Literal["ncrna", "cdna", "cds", "pep", "dna", "dna_index"],
    is_plant: bool,
) -> str: ...
def make_fasta_url(
    ensembl_release: Union[str, int],
    species: Union[str, Species],
    sequence_type: Literal["ncrna", "cdna", "cds", "pep", "dna", "dna_index"],
    is_plant: bool,
    server: str = ENSEMBL_FTP_SERVER,
    fasta_subdir=FASTA_SUBDIR_TEMPLATE,
) -> str: ...
