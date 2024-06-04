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

from typing import TYPE_CHECKING, List, Literal, Optional, Tuple, Union

from serializable import Serializable

if TYPE_CHECKING:
    from .database import Database
    from .exon import Exon
    from .gene import Gene
    from .locus import Locus
    from .sequence_data import SequenceData
    from .transcript import Transcript

class Genome(Serializable):
    """
    Bundles together the genomic annotation and sequence data associated with
    a particular genomic database source (e.g. a single Ensembl release) and
    provides a wide variety of helper methods for accessing this data.
    """

    def __init__(
        self,
        reference_name: str,
        annotation_name: str,
        annotation_version: Optional[Union[int, str]] = None,
        gtf_path_or_url: Optional[str] = None,
        transcript_fasta_paths_or_urls: Optional[list[str]] = None,
        protein_fasta_paths_or_urls: Optional[list[str]] = None,
        decompress_on_download: bool = False,
        copy_local_files_to_cache: bool = False,
        cache_directory_path: Optional[str] = None,
    ) -> None: ...
    @property
    def requires_gtf(self) -> bool: ...
    @property
    def requires_transcript_fasta(self) -> bool: ...
    @property
    def requires_protein_fasta(self) -> bool: ...
    def to_dict(self) -> dict: ...
    def _init_lazy_fields(self) -> None: ...
    def _get_cached_path(
        self, field_name, path_or_url, download_if_missing=False, overwrite=False
    ): ...
    def _get_gtf_path(self, download_if_missing=False, overwrite=False): ...
    def _get_transcript_fasta_paths(
        self, download_if_missing=False, overwrite=False
    ): ...
    def _get_protein_fasta_paths(self, download_if_missing=False, overwrite=False): ...
    def _set_local_paths(self, download_if_missing=True, overwrite=False): ...
    def required_local_files(self) -> List[str]: ...
    def required_local_files_exist(self, empty_files_ok: bool = False) -> bool: ...
    def download(self, overwrite: bool = False) -> None: ...
    def index(self, overwrite: bool = False) -> None: ...
    @property
    def db(self) -> Database: ...
    @property
    def protein_sequences(self) -> SequenceData: ...
    @property
    def transcript_sequences(self) -> SequenceData: ...
    def install_string(self) -> str: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
    def _fields(
        self,
    ) -> Tuple[str, str, Union[int, str], str, Tuple[str], Tuple[str]]: ...
    def __eq__(self, other) -> bool: ...
    def __hash__(self) -> int: ...
    def clear_cache(self) -> None: ...
    def delete_index_files(self) -> None: ...
    def _all_feature_values(
        self,
        column: str,
        feature: str,
        distinct: bool = True,
        contig: str = None,
        strand: Literal["+", "-"] = None,
    ) -> List: ...
    def transcript_sequence(self, transcript_id: str) -> Optional[str]: ...
    def protein_sequence(self, protein_id) -> Optional[str]: ...
    def genes_at_locus(
        self,
        contig: str,
        position: int,
        end: Optional[int] = None,
        strand: Optional[Literal["+", "-"]] = None,
    ) -> List[str]: ...
    def transcripts_at_locus(
        self,
        contig: str,
        position: int,
        end: Optional[int] = None,
        strand: Optional[Literal["+", "-"]] = None,
    ) -> List[str]: ...
    def exons_at_locus(
        self,
        contig: str,
        position: int,
        end: Optional[int] = None,
        strand: Optional[Literal["+", "-"]] = None,
    ) -> List[str]: ...
    def gene_ids_at_locus(
        self,
        contig: str,
        position: int,
        end: Optional[int] = None,
        strand: Optional[Literal["+", "-"]] = None,
    ) -> List[str]: ...
    def gene_names_at_locus(
        self,
        contig: str,
        position: int,
        end: Optional[int] = None,
        strand: Optional[Literal["+", "-"]] = None,
    ) -> List[str]: ...
    def exon_ids_at_locus(
        self,
        contig: str,
        position: int,
        end: Optional[int] = None,
        strand: Optional[Literal["+", "-"]] = None,
    ) -> List[str]: ...
    def transcript_ids_at_locus(
        self,
        contig: str,
        position: int,
        end: Optional[int] = None,
        strand: Optional[Literal["+", "-"]] = None,
    ) -> List[str]: ...
    def transcript_names_at_locus(
        self,
        contig: str,
        position: int,
        end: Optional[int] = None,
        strand: Optional[Literal["+", "-"]] = None,
    ) -> List[str]: ...
    def protein_ids_at_locus(
        self,
        contig: str,
        position: int,
        end: Optional[int] = None,
        strand: Optional[Literal["+", "-"]] = None,
    ) -> List[str]: ...

    ###################################################
    #
    #         Methods which return Locus objects
    #         containing (contig, start, stop, strand)
    #         of various genomic entities
    #
    ###################################################

    def locus_of_gene_id(self, gene_id: str) -> Locus: ...
    def loci_of_gene_names(self, gene_name: str) -> Locus: ...
    def locus_of_transcript_id(self, transcript_id: str) -> Locus: ...
    def locus_of_exon_id(self, exon_id: str) -> Locus: ...

    ###################################################
    #
    #                  Contigs
    #
    ###################################################

    def contigs(self) -> List[str]: ...

    ###################################################
    #
    #             Gene Info Objects
    #
    ###################################################

    def genes(
        self, contig: Optional[str] = None, strand: Optional[Literal["+", "-"]] = None
    ) -> List[Gene]: ...
    def gene_by_id(self, gene_id: str) -> Gene: ...
    def genes_by_name(self, gene_name: str) -> List[Gene]: ...
    def gene_by_protein_id(self, protein_id: str) -> Gene: ...

    ###################################################
    #
    #             Gene Names
    #
    ###################################################

    def _query_gene_name(
        self,
        property_name: str,
        property_value: str,
        feature_type: Literal["gene", "transcript", "exon"],
    ) -> str: ...
    def gene_names(
        self, contig: Optional[str] = None, strand: Optional[Literal["+", "-"]] = None
    ) -> List[str]: ...
    def gene_name_of_gene_id(self, gene_id: str) -> str: ...
    def gene_name_of_transcript_id(self, transcript_id: str) -> str: ...
    def gene_name_of_transcript_name(self, transcript_name: str) -> str: ...
    def gene_name_of_exon_id(self, exon_id: str) -> str: ...

    ###################################################
    #
    #             Gene IDs
    #
    ###################################################

    def _query_gene_ids(
        self,
        property_name: str,
        value: str,
        feature: Literal["gene", "CDS"] = "gene",
    ) -> List[str]: ...
    def gene_ids(
        self, contig: Optional[str] = None, strand: Optional[Literal["+", "-"]] = None
    ) -> List[str]: ...
    def gene_ids_of_gene_name(self, gene_name: str) -> str: ...
    def gene_id_of_protein_id(self, protein_id: str) -> str: ...

    ###################################################
    #
    #             Transcript Info Objects
    #
    ###################################################

    def transcripts(
        self, contig: Optional[str] = None, strand: Optional[Literal["+", "-"]] = None
    ) -> List[Transcript]: ...
    def transcript_by_id(self, transcript_id: str) -> Transcript: ...
    def transcripts_by_name(self, transcript_name: str) -> List[Transcript]: ...
    def transcript_by_protein_id(self, protein_id: str) -> Transcript: ...

    ###################################################
    #
    #            Transcript Names
    #
    ###################################################

    def _query_transcript_names(self, property_name: str, value: str) -> List[str]: ...
    def transcript_names(
        self, contig: Optional[str] = None, strand: Optional[Literal["+", "-"]] = None
    ) -> List[str]: ...
    def transcript_names_of_gene_name(self, gene_name: str) -> List[str]: ...
    def transcript_name_of_transcript_id(self, transcript_id: str) -> List[str]: ...

    ###################################################
    #
    #            Transcript IDs
    #
    ###################################################

    def _query_transcript_ids(
        self,
        property_name: str,
        value: str,
        feature: Literal["transcript", "CDS"] = "transcript",
    ) -> List[str]: ...
    def transcript_ids(
        self, contig: Optional[str] = None, strand: Optional[Literal["+", "-"]] = None
    ) -> List[str]: ...
    def transcript_ids_of_gene_id(self, gene_id: str) -> List[str]: ...
    def transcript_ids_of_gene_name(self, gene_name: str) -> List[str]: ...
    def transcript_ids_of_transcript_name(self, transcript_name: str) -> List[str]: ...
    def transcript_ids_of_exon_id(self, exon_id: str) -> List[str]: ...
    def transcript_id_of_protein_id(self, protein_id: str) -> List[str]: ...

    ###################################################
    #
    #             Exon Info Objects
    #
    ###################################################

    def exons(
        self, contig: Optional[str] = None, strand: Optional[Literal["+", "-"]] = None
    ) -> List[Exon]: ...
    def exon_by_id(self, exon_id: str) -> Exon: ...

    ###################################################
    #
    #                Exon IDs
    #
    ###################################################

    def _query_exon_ids(self, property_name: str, value: str) -> List[str]: ...
    def exon_ids(
        self, contig: Optional[str] = None, strand: Optional[Literal["+", "-"]] = None
    ) -> List[str]: ...
    def exon_ids_of_gene_id(self, gene_id: str) -> List[str]: ...
    def exon_ids_of_gene_name(self, gene_name: str) -> List[str]: ...
    def exon_ids_of_transcript_name(self, transcript_name: str) -> List[str]: ...
    def exon_ids_of_transcript_id(self, transcript_id: str) -> List[str]: ...
    def protein_ids(
        self, contig: Optional[str] = None, strand: Optional[Literal["+", "-"]] = None
    ) -> List[str]: ...
