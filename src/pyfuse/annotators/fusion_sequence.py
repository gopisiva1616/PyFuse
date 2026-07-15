#!/usr/bin/env python3

"""Retrieve fusion sequences from teh reference fastq file based on the input fusion coordinates"""

import os
import sys
import csv
import time
import pyfaidx
import logging
import argparse
import progressbar
import pandas as pd
from Bio.Seq import Seq
from Bio.Data import CodonTable
from pyfuse.utils.common_utils import utils, config

start_time = time.time()
logger = logging.getLogger(__name__)
csv.field_size_limit(sys.maxsize)


def validatefilepaths(path):
    """
    Function to check for validatiy of file paths provided

    Attrbutes
    ---------
    path to a certain directory of file : str

    """

    if os.path.exists(path) == False:
        raise Exception("ERROR: The given path" + path + "does not exist or the path is incorrect \n Please check the path")
    if os.path.getsize(path) == 0:
        raise Exception("ERROR: The given path" + path + "does not exist or is inaccessible")
    return path


def argumentparser():
    help_text = "Extract fusion sequence for a given fusion candidate"
    parser = argparse.argumentparser(description=help_text, epilog=help_text)
    parser.add_argument("-i", help="input fusion file with first column as fusion name \
                        example, EML4-ALK", required="True", type=ValidatefilePaths)
    parser.add_argument("-r", help="Reference fasta file")
    parser.add_argument("-o", help="Path to output directory", type=ValidatefilePaths)
    parser.add_argument("-c", help="Coordinate file with exon start stop")
    return parser


def bedtools_get_fasta(chrs, seq_start, seq_end, strand):
    global pyref
    pyseq = pyref[chrs][int(seq_start):int(seq_end)]

    if strand == "-":
        pyseq = pyseq.reverse.complement

    return str(pyseq)


def nt_to_peptide_old(nucleotide_seq: str, frame: int = 1, to_stop: bool = True) -> str:
    """
    Convert a nucleotide sequence to a peptide sequence.

    Parameters:
    - nucleotide_seq (str): The input nucleotide sequence (A, T, C, G).
    - frame (int): Reading frame (1, 2, or 3).
    - to_stop (bool): If True, translation will stop at the first stop codon.

    Returns:
    - str: Peptide (amino acid) sequence.
    """
    # Clean sequence
    nucleotide_seq = nucleotide_seq.upper().replace("\n", "").replace(" ", "").replace("U", "T").replace("-", "")

    if nucleotide_seq == "NA":
        return "NA"

    if frame not in [1, 2, 3]:
        raise ValueError("Frame must be 1, 2, or 3.")

    # Adjust sequence for reading frame
    #seq = Seq(nucleotide_seq[frame - 1:])

    try:
        print('\n\n----Nucleotide_seq:', nucleotide_seq)
        peptide = nucleotide_seq.translate(to_stop=True)
    except Exception as e:
        raise ValueError(f"Translation error: {e}")
    
    print('Peptide sequence:', peptide)

    return str(peptide)

from Bio.Seq import Seq

def nt_to_peptide(nucleotide_seq: str, frame: int = 1, to_stop: bool = True) -> str:
    """
    Convert a nucleotide sequence to a peptide sequence.

    Parameters:
    - nucleotide_seq (str): The input nucleotide sequence (A, T, C, G).
    - frame (int): Reading frame (1, 2, or 3).
    - to_stop (bool): If True, translation will stop at the first stop codon.

    Returns:
    - str: Peptide (amino acid) sequence.
    """
    # Clean sequence
    nucleotide_seq = nucleotide_seq.upper().replace("\n", "").replace(" ", "").replace("U", "T").replace("-", "")

    if nucleotide_seq == "NA":
        return "NA"

    if frame not in [1, 2, 3]:
        raise ValueError("Frame must be 1, 2, or 3.")

    # Create a Seq object and apply reading frame
    seq = Seq(nucleotide_seq[frame - 1:])

    try:
        #print('\n\n----Nucleotide_seq:', seq)
        peptide = seq.translate(to_stop=to_stop)
    except Exception as e:
        raise ValueError(f"Translation error: {e}")

    #print('Peptide sequence:', peptide)

    return str(peptide)


def fetch_details_of_previous_exon(gt_combo, exon):
    """
    Function to fetch details of the preceeding exon from the
    exon of the breakpoint

    Attributes
    ----------

    gt_combo : Key of gene:transcript (of the fusion candidate) consisting of values as
     (start,stop, exon number etc..)
    """
    for keys in dict_start_stop:
        if (gt_combo == keys):
            for index, item in enumerate(dict_start_stop[gt_combo]):
                if (item[4] == exon):
                    neigboring_exon_element = dict_start_stop[gt_combo][index-1]

    return neigboring_exon_element


def fetch_details_of_next_exon(gt_combo, exon):
    for keys in dict_start_stop:
        if (gt_combo == keys):
            for index, item in enumerate(dict_start_stop[gt_combo]):
                if (item[4] == exon):
                    neigboring_exon_element = dict_start_stop[gt_combo][index+1]

    return neigboring_exon_element


# Retrive sequence for 5 prime end of the gene
def get_seq_5p(gt_combo_5p, strand_5p, exon_5p, exon_loc_5p, chr_5p, coord_5p):
    """
    Function to get sequence of the 5 prime candidate

    Attributes
    ----------
    fus_candidate : list of variables consisting of meta data for a
    given fusion breakpoint such as gene, transcript, exon , strand.
    """

    fasta = 'NA'

    if gt_combo_5p in dict_start_stop:
        for index, j in enumerate(dict_start_stop[gt_combo_5p]):
            if (str(j[4]) == exon_5p):
                chrs, seq_start, seq_end = j[0], int(j[1]), int(j[2])
                break  

        if exon_loc_5p == "Body":
            if strand_5p == '+':
                seq_start = coord_5p - 200
                seq_end = coord_5p                  
            elif strand_5p == '-':
                seq_start = coord_5p
                seq_end = coord_5p + 200

        if exon_loc_5p == "Start":
            if len(dict_start_stop[gt_combo_5p]) == 1:
                return 'NA'
            neigboring_exon_element = fetch_details_of_previous_exon(gt_combo_5p,
                                                                        exon_5p)
            seq_start = int(neigboring_exon_element[1])-1
            seq_end = neigboring_exon_element[2]

        if exon_loc_5p == "End":
            seq_start = seq_start-1

        fasta = bedtools_get_fasta(chrs, seq_start, seq_end, strand_5p)
        return fasta

# Retrieve sequence for 3' prime end of the gene
def get_seq_3p(gt_combo_3p, strand_3p, exon_3p, exon_loc_3p, chr_3p, coord_3p):
    """
    Function to get sequence of the 3 prime candidate

    Attributes
    ----------
    fus_candidate : list of variables consisting of meta data for a
    given fusion breakpoint such as gene, transcript, exon , strand.
    """

    fasta = 'NA'
    gene = gt_combo_3p.split(":")[0]
    transcript = gt_combo_3p.split(":")[1]

    if gt_combo_3p in dict_start_stop:
        for index, j in enumerate(dict_start_stop[gt_combo_3p]):
            if (str(j[4]) == exon_3p):
                chrs, seq_start, seq_end = j[0], int(j[1]), int(j[2])
                break
            
        if exon_loc_3p == "Body":
            if strand_3p == '+':
                seq_start = coord_3p
                seq_end = coord_3p + 200
            elif strand_3p == '-':        
                seq_start = coord_3p - 200
                seq_end = coord_3p                
                
        if exon_loc_3p == "End":
            if len(dict_start_stop[gt_combo_3p]) == 1:
                return 'NA'
            neigboring_exon_element = fetch_details_of_next_exon(gt_combo_3p,
                                                                    exon_3p)
            seq_start = int(neigboring_exon_element[1])-1
            seq_end = neigboring_exon_element[2]
        
        if exon_loc_3p == "Start":
            seq_start = seq_start-1
        
        fasta = bedtools_get_fasta(chrs, seq_start, seq_end, strand_3p)
        return fasta


def parse_coodinates(j, header_index):
    """
    Function to parse information on 5 prime and 3 prime coordinates
    """
    gene = j[header_index].split("|")[2]
    transcript = j[header_index].split("|")[3]
    gt_combo = str(gene) + ":" + str(transcript)
    strand = j[header_index].split("|")[0]
    exon = str(j[header_index].split("|")[1].split("_")[1])
    exon_loc = str(j[header_index].split("|")[1].split("_")[0])
    return gt_combo, strand, exon, exon_loc


def main():
    parser = argumentparser()
    try:
        args = parser.parse_args()
    except:
        sys.exc_info()[1]
        print(f'Please provide above given arguments')
        exit()

    input_file = args.i


def fusion_sequence(input_fusions, ref_file):
    """
    Function to retrieve fusion seqeunce by retrieving sequence of the
    exons of teh 5 prime and 3 prime gene partners

    Attributes
    ----------
    reference - Ref.fa for hg19 for extraction fasta seqeunces of the fusion partners.
    fusions_file - Input data frame coming from fusion annotation code
    exon_file - Reference bed file consisting of exon start stop, exon number and strand
    information
    """

    global reference, pyref
    logger.info("-- Determining the fusion sequence")


    fusions_file = utils.read_df_or_file(input_fusions)
    reference = ref_file
    pyref = pyfaidx.Fasta(reference)
    exon_file = utils.read_df_or_file(config['all_gene_exons'])

    # From the exon list file make a dictionary of gene transcript coordinates
    global dict_start_stop
    dict_start_stop = {}

    for count, line in enumerate(exon_file):
        if count == 0:
            continue
        i = line.strip().split('\t')
        my_gene = i[4].split("|")[2]
        my_transcript = i[4].split("|")[3]
        my_key = my_gene + ":" + my_transcript
        direction = i[4].split("|")[0]
        ex = i[4].split("|")[1]
        my_value = [i[0], i[1], i[2], direction, ex]

        if my_key not in dict_start_stop:
            dict_start_stop[my_key] = [my_value]
        else:
            dict_start_stop[my_key].append(my_value)

    final_list = []
    bar = progressbar.ProgressBar(max_value=len(input_fusions))
    for idx, line in enumerate(fusions_file):

        try:

            final_seq = 'NA'
            line = line.replace('"', '')
            j = str(line).strip().split('\t')

            if idx == 0:
                header = j
                header.append("Fusion_nucleotide_sequence")
                header.append("Fusion_peptide_sequence")
                #print('\nFUS header list', header, type(header))

                header_5p_index = j.index("5'_Exon_Annotation")
                coordinate_5p = j.index("5'co-ordinate")
                header_3p_index = j.index("3'_Exon_Annotation")
                coordinate_3p = j.index("3'co-ordinate")
                continue
                # 5' prime exon annotation
            else:

                if (j[header_5p_index] != "NA"):
                    gt_combo_5p, strand_5p, exon_5p, exon_loc_5p = parse_coodinates(j, header_5p_index)
                    chr_5p = j[coordinate_5p].split(":")[0]
                    coord_5p = int(j[coordinate_5p].split(":")[1])
                    seq_5p = get_seq_5p(gt_combo_5p, strand_5p, exon_5p, exon_loc_5p, chr_5p, coord_5p)

                else:

                    chrs_5p = j[coordinate_5p].split(":")[0]
                    coord_5p = int(j[coordinate_5p].split(":")[1])
                    seq_start_5p = coord_5p - 200
                    seq_end_5p = coord_5p
                    strand_5p = "NA"
                    seq_5p = bedtools_get_fasta(chrs_5p, seq_start_5p, seq_end_5p, strand_5p)

                # 3' prime exon annotation
                if (str(j[header_3p_index]) != "NA"):
                    gt_combo_3p, strand_3p, exon_3p, exon_loc_3p = parse_coodinates(j, header_3p_index)
                    chr_3p = j[coordinate_3p].split(":")[0]
                    coord_3p = int(j[coordinate_3p].split(":")[1])
                    seq_3p = get_seq_3p(gt_combo_3p, strand_3p, exon_3p,
                                        exon_loc_3p, chr_3p, coord_3p)
                else:
                    #print('\n\n----3 prime exon annotation:', j[header_3p_index])


                    chrs_3p = j[coordinate_3p].split(":")[0]
                    coord_3p = int(j[coordinate_3p].split(":")[1])
                    seq_start_3p = coord_3p - 1
                    seq_end_3p = coord_3p + 200
                    strand_3p = "NA"
                    seq_3p = bedtools_get_fasta(chrs_3p, seq_start_3p, seq_end_3p, strand_3p)
            try:
                final_seq = seq_5p + "-" + seq_3p
            except TypeError as e:
                logger.error(f'\nTYPE-ERROR:  {e} out,{seq_5p} - {seq_3p}  \n{line}')

        except IndexError:
            logger.error(f'\nFAILING INPUT lINE: {line}')

        
        peptide_seq = nt_to_peptide(final_seq)
        j.append(final_seq)
        j.append(peptide_seq)
        final_list.append(j)
        bar.update(idx)

    final_df = pd.DataFrame(final_list)
    final_df.columns = header
    #final_df.to_csv('finaldf.csv', header=True, sep='\t')
    return final_df


if __name__ == '__main__':
    main()
