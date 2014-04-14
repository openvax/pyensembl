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
from Bio.Seq import Seq

from annotation_data import EnsemblAnnotationData

data = EnsemblAnnotationData()

def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())

def get_exons_from_transcript(transcript_id):
    """
    Filter exons down to those with this transcript_id

    Parameters
    ----------
    transcript_id : transcript id, of the from EST#####

    returns Pandas dataframe containing only those exons
    """
    exon_data = data.exon_data
    exons = exon_data[exon_data['stable_id_transcript'] == transcript_id]
    assert len(exons) > 0, \
        "Couldn't find exons for transcript %s" % transcript_id
    fields = [
        'exon_id',
        'seq_start',
        'start_exon_id',
        'seq_end',
        'end_exon_id',
        'stable_id_exon',
        'seq_region_start_exon',
        'seq_region_end_exon'
    ]
    return exons[fields]


def get_idx_from_interval(pos, intervals):
    idx = 0
    for (start, end) in intervals:
        if pos > end:
            idx += (end - start) + 1
        elif pos <= end and pos >= start:
            return idx + (pos - start)
        else:
            return None

def get_strand(transcript_id):
    """
    Gets the strand of the gene

    Parameters
    ----------
    transcript_id :
        Transcript id, of the from EST#####

    Return strand : int, +1 for forward strand else -1
    """
    transcript_data = data.transcript_data
    transcript_info = transcript_data[transcript_data['stable_id_transcript'] == transcript_id]

    if transcript_info.empty:
        logging.warn("Transcript %s has no sequence information", transcript_id)
        return 1

    strand = list(transcript_info['seq_region_strand_gene'])[0]
    return strand

def is_forward_strand(transcript_id):
    return get_strand(transcript_id) > 0

def is_incomplete_cds(transcript_id):
    """
    Compute 5 prime incomplete CDS - checks the start_phase of the first exon

    transcript_id :
        Transcript id, of the from EST#####
    """
    start_phase = get_cds_start_phase(transcript_id)
    if start_phase:
        return start_phase > 0
    else:
        return False

def get_start_exon(transcript_id):
    exon_data = data.exon_data
    first_exon_mask = (exon_data['stable_id_transcript'] == transcript_id) & (exon_data['rank'] == 1)
    start_exon = exon_data[first_exon_mask]
    if start_exon.empty:
        return None
    else:
        return start_exon.to_dict(outtype = 'records')[0]

def get_cds_start_phase(transcript_id):
    """
    Compute CDS start_phase - checks the phase of the first exon

    transcript_id :
        Transcript id, of the from EST#####
    """
    start_exon = get_start_exon(transcript_id)
    if start_exon is not None:
        return max(start_exon['phase'], 0)  # -1 means no phase (equiv. to zero)
    else:
        return None

def get_transcript_index_from_pos(
        pos, 
        transcript_id,
        skip_untranslated_region= True):
    """
    Gets the index into to the transcript from genomic position
    The transcript is composed of spliced exons that have genomic start and
    stop positions.  This function searches for the exon that matches this
    genomic positions and returns the index into the transcript

    Parameters
    ----------
    position : int
        Genomic position in the contig

    transcript_id :
        Transcript id, of the from EST#####

    skip_untranslated_region : bool, optional
        If True (default), then give position in the CDS (coding sequence),
        otherwise give position in the longer full cDNA sequence.
    """
    exons = get_exons_from_transcript(transcript_id)
    exons = exons.sort(columns=['seq_region_start_exon', 'seq_region_end_exon'])
    exons['exon_length'] = exons['seq_region_end_exon'] - exons['seq_region_start_exon'] + 1
    starts = exons['seq_region_start_exon']
    stops = exons['seq_region_end_exon']
    intervals = zip(starts, stops)

    transcript_length = exons['exon_length'].sum()
    logging.info("Full transcript length for %s = %d", transcript_id, transcript_length)
    transcript_idx = get_idx_from_interval(pos, intervals)

    if transcript_idx is None:
        logging.warning("Couldn't find position %d in transcript %s",
            pos, transcript_id)
    else:
        # Reverse array index if on reverse strand
        forward = is_forward_strand(transcript_id)
        transcript_idx = transcript_idx if forward else transcript_length - transcript_idx - 1
        logging.info("Transcript strand forward? %s = %d", transcript_id, forward)
        if skip_untranslated_region:
            # Adjust for translations (CDS) start region
            prefix_utr_length = get_five_prime_utr_length(exons, forward)
            logging.info("UTR length for %s = %d", transcript_id, prefix_utr_length)
            if transcript_idx < prefix_utr_length:
                logging.warn("UTR mutation at cDNA position %d, transcript %s does not apply", 
                    transcript_idx, transcript_id)
                return None
            else:
                transcript_idx -= prefix_utr_length

        # Adjust for CDS start phase if first exon is out of phase
        transcript_phase = get_cds_start_phase(transcript_id)
        transcript_idx += transcript_phase
        if transcript_phase > 0:
            logging.warn("Transcript %s is incomplete", transcript_id)

        # TODO: check that index is within the mRNA transcript
        # need to get the length of the coding region from the transcript_id
        #suffix_utr_length = get_three_prime_utr_length(exons, forward)
        #assert transcript_idx <= transcript_length + suffix_utr_length

    return transcript_idx

def get_five_prime_utr_length(exons_df, forward = True):
    """
    Gets the length of the 5' UTR from a set of sorted exons from a specifc transcript

    Parameters
    ----------
    exons : Pandas dataframe with 'exon_id', 'seq_region_end_exon', 'seq_region_start_exon'
            Also, 'start_exon_id' marks the exon that starts translation and 'seq_start'
            is the offset into the first translated exon

    forward : bool, default = True, is forward strand or not

    Return utr_length : int
    """
    exons_df = exons_df.sort(columns=['seq_region_start_exon', 'seq_region_end_exon'], ascending=[forward, forward])
    utr_length = 0
    for idx, row in exons_df.iterrows():
        if row['exon_id'] == row['start_exon_id']:
            utr_length += row['seq_start'] - 1
            return utr_length         
        else:
            utr_length += row['seq_region_end_exon'] - row['seq_region_start_exon'] + 1
    return None

def get_three_prime_utr_length(exons_df, forward = True):
    """
    Gets the length of the 3' UTR from a set of sorted exons from a specifc transcript

    Parameters
    ----------
    exons : Pandas dataframe with 'exon_id', 'seq_region_end_exon', 'seq_region_start_exon'
            Also, 'end_exon_id' marks the exon that starts translation and 'seq_end'
            is the offset into the last translated exon

    forward : bool, default = True, is forward strand or not

    Return utr_length : int
    """
    reverse = not forward
    exons_df = exons_df.sort(columns=['seq_region_start_exon', 'seq_region_end_exon'], ascending=[reverse, reverse])
    utr_length = 0
    for idx, row in exons_df.iterrows():
        exon_length = row['seq_region_end_exon'] - row['seq_region_start_exon'] + 1
        if row['exon_id'] == row['end_exon_id']:
            utr_length += exon_length - row['seq_end']
            return utr_length           
        else:
            utr_length += exon_length
    return None

def annotate(
        vcf_df,
        annotation_df,
        predicate,
        left_col = 'chr',
        right_col = 'name'):
    crossed = vcf_df.merge(
        annotation_df, left_on=left_col, right_on=right_col, how='left')
    annotated = crossed[crossed.apply(predicate, axis=1)]
    return annotated.drop_duplicates()


def _transcript_matches(transcript_row):
    position = transcript_row['pos']
    contig = transcript_row['chr']
    return transcript_row['seq_region_start_transcript'] < position \
            and transcript_row['seq_region_end_transcript'] > position \
            and transcript_row['name'] == contig

def annotate_vcf_transcripts(vcf_df):
    """
    Get list of transcript id from position

    Parameters
    ----------
    vcf : Pandas dataframe with chr, pos, ref, alt columns

    Return df with extra columns:
        'name', 'stable_id_gene', 'description_gene',
        'seq_region_start_gene', 'seq_region_end_gene',
        'stable_id_transcript', 'seq_region_start_transcript',
        'seq_region_end_transcript'

    """
    if 'gene_stable_id' in vcf_df.columns:
        annotated = annotate(vcf_df,
            data.transcript_data,
            _transcript_matches,
            left_col=['chr', 'gene_stable_id'],
            right=['name', 'gene_stable_id'])
    else:
        annotated = annotate(vcf_df, data.transcript_data, _transcript_matches)
    return annotated



def _gene_matches(gene_row):
    position = gene_row['pos']
    contig = gene_row['chr']
    return gene_row['seq_region_start_gene'] <= position \
            and gene_row['seq_region_end_gene'] > position \
            and gene_row['name'] == contig


def annotate_vcf_genes(vcf_df):
    """
    Get list of gene id from position

    Parameters
    ----------
    vcf_df : Pandas dataframe
        Expected to have columns 'chr', 'pos', 'ref', 'alt'

    Return df with extra columns:
        'name', 'stable_id_gene', 'description_gene',
        'seq_region_start_gene', 'seq_region_end_gene'

    """
    return annotate(vcf_df, data.gene_data, _gene_matches)

