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

from typing import TYPE_CHECKING, Any, List, Literal, Optional, Union
from .common import memoize

if TYPE_CHECKING:
    import logging
    import polars
    from .locus import Locus
    from sqlite3 import Connection

# any time we update the database schema, increment this version number
DATABASE_SCHEMA_VERSION = 3

logger: logging.Logger = ...

class Database(object):
    def __init__(
        self,
        gtf_path: str,
        install_string: Optional[str] = None,
        cache_directory_path: Optional[str] = None,
        restrict_gtf_columns: Optional[List[str]] = None,
        restrict_gtf_features: Optional[List[str]] = None,
    ) -> None: ...
    def __eq__(self, other) -> bool: ...
    def __str__(self) -> str: ...
    def __hash__(self) -> int: ...
    @property
    def local_db_filename(self) -> str: ...
    @property
    def local_db_path(self) -> str: ...
    def _all_possible_indices(self, column_names: str) -> List[List[str]]: ...

    PRIMARY_KEY_COLUMNS = {"gene": "gene_id", "transcript": "transcript_id"}

    def _get_primary_key(
        self, feature_name: str, feature_df: polars.DataFrame
    ) -> str: ...
    def _feature_indices(
        self, all_index_groups: List, primary_key: str, feature_df: polars.DataFrame
    ) -> List: ...
    def create(self, overwrite: bool = False) -> Connection: ...
    def _get_connection(self) -> Connection: ...
    @property
    def connection(self) -> Connection: ...
    def connect_or_create(self, overwrite: bool = False) -> Connection: ...
    def columns(self, table_name: str) -> List[str]: ...
    def column_exists(self, table_name: str, column_name: str) -> bool: ...
    def column_values_at_locus(
        self,
        column_name: str,
        feature: str,
        contig: str,
        position: int,
        end: Optional[int] = None,
        strand: Literal["+", "-"] = None,
        distinct: bool = False,
        sorted: bool = False,
    ) -> List[Any]: ...
    def distinct_column_values_at_locus(
        self,
        column: str,
        feature: str,
        contig: str,
        position: int,
        end: Optional[int] = None,
        strand: Literal["+", "-"] = None,
    ) -> List[Any]: ...
    def run_sql_query(
        self, sql: str, required: bool = False, query_params: List[Union[str, int]] = []
    ) -> List[Any]: ...
    @memoize
    def query(
        self,
        select_column_names: List[str],
        filter_column: str,
        filter_value: str,
        feature: Literal["transcript", "gene", "exon", "CDS"],
        distinct: bool = False,
        required: bool = False,
    ) -> List[Any]: ...
    def query_one(
        self,
        select_column_names: List[str],
        filter_column: str,
        filter_value: str,
        feature: Literal["transcript", "gene", "exon", "CDS"],
        distinct: bool = False,
        required: bool = False,
    ): ...
    @memoize
    def query_feature_values(
        self,
        column: str,
        feature: Literal["transcript", "gene", "exon", "CDS"],
        distinct: bool = True,
        contig: Optional[str] = None,
        strand: Optional[Literal["+", "-"]] = None,
    ) -> List[str]: ...
    def query_distinct_on_contig(
        self,
        column_name: str,
        feature: Literal["transcript", "gene", "exon", "CDS"],
        contig: str,
    ) -> List[str]: ...
    def query_loci(
        self,
        filter_column: str,
        filter_value: str,
        feature: Literal["transcript", "gene", "exon", "CDS"],
    ) -> List[Locus]: ...
    def query_locus(
        self,
        filter_column: str,
        filter_value: str,
        feature: Literal["transcript", "gene", "exon", "CDS"],
    ) -> Locus: ...
    def _load_gtf_as_dataframe(
        self, usecols: Optional[List[str]] = None, features: Optional[List[str]] = None
    ) -> polars.DataFrame: ...
