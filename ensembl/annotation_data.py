# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging 

import pandas as pd

from transcript_metadata import download_transcript_metadata

class EnsemblAnnotationData(object):

    """
    Singleton class which allows for lazy loading of
    exon/transcript annotations
    """

    def __init__(self):
        pass

    @property
    def transcript_metadata_path(self):
        if not hasattr(self, '_transcript_metadata_path'):
            self._transcript_metadata_path = download_transcript_metadata()
        return self._transcript_metadata_path

    @property
    def exon_data(self):
        """
        Dataframe containing exon data
        """
        if not hasattr(self, '_exon_data'):
            path = self.transcript_metadata_path
            self._exon_data = pd.read_csv(path, sep='\t')
            logging.info("Loaded exon metadata with columns %s", list(self._exon_data.columns)) 
        return self._exon_data

    @property
    def transcript_data(self):
        """
        Subset columns for transcript data only
        """
        if not hasattr(self, '_transcript_data'):
            transcript_cols = [
                'name', 'stable_id_gene', 'description_gene',
                'seq_region_start_gene', 'seq_region_end_gene',
                'seq_region_strand_gene', 'stable_id_transcript', 
                'seq_region_start_transcript', 'seq_region_end_transcript'
            ]
            self._transcript_data = self.exon_data[transcript_cols].drop_duplicates()
        return self._transcript_data

    @property
    def gene_data(self):
        """
        Subset columns for gene data only
        """
        if not hasattr(self, '_gene_data'):
            gene_cols = [
                'name', 'stable_id_gene', 'description_gene',
                'seq_region_start_gene', 'seq_region_end_gene'
            ]
            self._gene_data = self.transcript_data[gene_cols].drop_duplicates()
        return self._gene_data
