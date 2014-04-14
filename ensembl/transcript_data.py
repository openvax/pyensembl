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
from os.path import join
import sqlite3

from Bio import SeqIO
import appdirs
from epitopes.download import fetch_fasta_db, ensure_dir

CDNA_TRANSCRIPT_URL = \
'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz'

CDNA_TRANSCRIPT_FILE = 'Homo_sapiens.GRCh37.75.cdna.all.fa'

def _build_cdna_db():
    """
    Download a FASTA file containing cDNA sequences for each known transcript,
    return sqlite3 database mapping ensembl transcript IDs to cDNA sequences.
    """
    return fetch_fasta_db(
        table_name = "CDNA",
        fasta_filename = CDNA_TRANSCRIPT_FILE,
        download_url = CDNA_TRANSCRIPT_URL,
        key_column = 'id',
        value_column = 'seq',
        subdir = "immuno")

CDS_TRANSCRIPT_URL = \
'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.75.cds.all.fa.gz'

CDS_TRANSCRIPT_FILE = 'Homo_sapiens.GRCh37.75.cds.all.fa'

def _build_cds_db():
    """
    Download a FASTA file containing CDS sequences for each known transcript,
    return sqlite3 database mapping ensembl transcript IDs to CDS sequences.
    """
    return fetch_fasta_db(
        table_name = "CDS",
        fasta_filename = CDS_TRANSCRIPT_FILE,
        download_url = CDS_TRANSCRIPT_URL,
        key_column = 'id',
        value_column = 'seq',
        subdir = "immuno")

PROTEIN_TRANSCIPT_URL = \
'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.75.pep.all.fa.gz'

PROTEIN_TRANSCRIPT_FILE = 'Homo_sapiens.GRCh37.75.pep.all.fa'

def _build_protein_db():
    """
    Downloads a FASTA file containing amino acid sequences of human
    transcripts, return sqlite3 database mapping ensembl protein IDs to
    amino acid sequences.
    """
    return fetch_fasta_db(
        table_name = "PROTEIN",
        fasta_filename = PROTEIN_TRANSCRIPT_FILE,
        download_url = PROTEIN_TRANSCIPT_URL,
        key_column = 'id',
        value_column = 'seq',
        subdir = "immuno")

def _exec_transcript_query(db, table_name, transcript_id):
    query = "select seq from %s where id = ?" % table_name

    cursor = db.execute(query, (transcript_id,))
    results = cursor.fetchmany()
    if len(results) == 0:
        logging.warning("No entries found with transcript_id = %s",
            transcript_id)
        return None
    else:
        assert len(results) == 1, \
            "Too many entries (%d) with transcript_id = %s" % \
            (len(results), transcript_id)
        # get the first result
        # and return the first element of its tuple
        return results[0][0]

class EnsemblReferenceData(object):
    """
    Singleton class which allows for lazy loading of reference
    cDNA and amino acid sequences of transcripts
    """
    def __init__(self):
        self._cdna_db = None
        self._cds_db = None
        self._protein_db = None

    def get_cdna(self, transcript_id):
        if self._cdna_db is None:
            self._cdna_db = _build_cdna_db()
        return _exec_transcript_query(
            self._cdna_db, "CDNA", transcript_id)

    def get_cds(self, transcript_id):
        if self._cds_db is None:
            self._cds_db = _build_cds_db()
        return _exec_transcript_query(
            self._cds_db, "CDS", transcript_id)

    def get_protein(self, transcript_id):
        if self._protein_db is None:
            self._protein_db = _build_protein_db()
        print transcript_id
        return _exec_transcript_query(
            self._protein_db, "PROTEIN", transcript_id)
