#!/usr/bin/env python3

"""Calculate frame of a given fusion candidate"""

import os
import sys
import csv
import time
import logging
import argparse
import pandas as pd
from tqdm import tqdm
from pyfuse.utils.common_utils import utils, config

start_time = time.time()
logger = logging.getLogger(__name__)
csv.field_size_limit(sys.maxsize)


def ArgumentParser():
    """
    Function to parse arguments from the command line
    Arguments include
    1. Input to a file listing fusion candidates
    2. Path to output directory
    3. Path to tab delimited text file consisting of information on
       chr , start codon positio, start of exon, stop of exon, exon length,
       exon number, gene and transcript.
    """

    help_text = """Calculates the frame of an incoming fusion candidate"""
    parser = argparse.ArgumentParser(description=help_text, epilog=help_text)
    parser.add_argument("-i", help="Path to incoming fusion file, EML4-ALK",
                        required="True", type=ValidatefilePaths)
    parser.add_argument("-o", help="Path to output directory", default=os.getcwd())
    parser.add_argument("-c", help="Path to start codon file")
    return parser


def get_frame_status(frame_5p, frame_3p):
    """
    Function to calculate the frame status of a given fusion candidate using
    the values of the five prime end of the fusion and the three prime end of the
    fusion.

    The function will provide a set of conditions for the various combinations
    of 5' and 3 ' frame status values being output from the compute_frame
    function and provide a value for the frame status of the fusion
    candidate using the respective combinations.

    Attributes
    ----------
    frame_5p, frame_3p: str
        Frame values for 5' and 3' partners of a fusion candidate.
        Acceptable frame values for both 5' and 3' gene partners include:
        1 , 2, 3, In_Untranslated_Region, Start_Codon_NA, NA
    """
    frame_5p = str(frame_5p)
    frame_3p = str(frame_3p)

    utr_scd_na = ["In_Untranslated_Region", "Start_Codon_NA"]
    in_frame_combos = [("1", "2"), ("2", "3"), ("3", "1")]
    out_frame_combos = [("1", "1"), ("1", "3"), ("2", "1"),
                        ("2", "2"), ("3", "2"), ("3", "3")
                        ]

    if (frame_5p, frame_3p) in in_frame_combos:
        frame_status = "In-Frame"

    elif (frame_5p, frame_3p) in out_frame_combos:
        frame_status = "Out-of-Frame"

    elif frame_5p in ["In_Untranslated_Region", "Start_Codon_NA", 'NA'] and \
            frame_3p in ['1', '2', '3']:
        frame_status = "5UTR-CDS"

    elif frame_3p in ["In_Untranslated_Region", "Start_Codon_NA", 'NA'] and \
            frame_5p in ['1', '2', '3']:
        frame_status = "CDS-3UTR"

    elif frame_5p in utr_scd_na and frame_3p in utr_scd_na:
        frame_status = "5UTR-3UTR"

    elif frame_3p in utr_scd_na and frame_5p == 'NA':
        frame_status = "INTRON-3UTR"

    elif frame_5p in utr_scd_na and frame_3p == 'NA':
        frame_status = "5UTR-INTRON"

    elif frame_5p == 'NA' and frame_3p == 'NA':
        frame_status = "INTRON-INTRON"
    else:
        frame_status = "UNK"

    return frame_status

results_list = []
def compute_frame(exon_sent, exon_name, extra_bases, gene_trans,
                  index_of_bp_exon, bkpt_coordinate, condition):
    """
    Function to compute the frame of a gene partner involved in a fusion.
    This function is called for both 5' and 3' partners of a fusion ;
    first for 5' prime and then for 3' prime.

    1. The function first tries to locate the position of the start codon (SC)
    within a given transcript; i.e it tries to find if the start codon is
    present within the first exon or any of the other exons of a transcript.

    2. Post that the function tries to then locate the position of the start
    codon within the exon consisting of the start codon and uses different
    conditions for determining the coding length of the exon to be included
    in the calculation for the frame length using both negative and positive strands

    Attributes
    ----------
    exon_sent, gene_trans: str
        1. Exon consisting of the breakpoint
        2. gene_trans : A key consisting of the gene:transcript of the gene partner

     exon_name, extra_bases, exon_len, index_of_bp_exon : int
        1. exon_name is the index used to iterate through all exons of a particular transcript.
           Each exon_name is composed of chromosome, start codon position, exon start, exon end
           length of an exon, strand
        2. extra_bases is the number of bases to be added if the breakpoint lies
           within the exon
        3. exon_len is the length of the exon consisting of the breakpoint
        4. index_of_bp_exon is the index of the exon consisting of the breakpoint.

     The function returns the value of the frame of the 5' gene partner or the
     3' gene partner.
    """
    final_len, length_intermed_exons, modulo, frame = 0, 0, 0, 0
    result_dict = {
        'gene_transcript': f'{gene_trans}:{exon_sent}',
        'coord': bkpt_coordinate,
        'condition': condition,
        'length_intermediate': length_intermed_exons, 
        'extraBases': extra_bases,
        'final_length': final_len,
        'modulo': modulo,
        'frame': frame
    }
    bp_exon = int(exon_sent.split("E")[1]) 
    strand = str(start_codon_dict[gene_trans][0][5])


    # if the position of the SC in the start codon file is not "-"
    if exon_name[1] != "-":
        start_codon = int(exon_name[1])

    if start_codon_dict[gene_trans][index_of_bp_exon][1] == "-":

        return result_dict, "In_Untranslated_Region"
    if start_codon_dict[gene_trans][index_of_bp_exon][1] == "0":
        return result_dict, "Start_Codon_NA"
    else:
        # Find the location of the SC in the list of exons for a transcript


        start_codon_list = []
        my_list = start_codon_dict[gene_trans]
        for a1 in range(0, len(my_list)):
            temp = my_list[a1][1]
            start_codon_list.append(temp)
        start_codon_idx = start_codon_list.index(str(start_codon)) #start_codon_exon
    
    # var assignmnets
    start_codon_exon_start, start_codon_exon_end = 0, 0
    bkpt_exn_len = abs(extra_bases)
    start_codon_pos_list = start_codon_dict[gene_trans][start_codon_idx]
    start_codon_exn_len = int(start_codon_dict[gene_trans][start_codon_idx][4])
    start_codon_to_exn_end = abs(start_codon - int(start_codon_exon_end)) + 1  # due to subtraction, adding 1
    start_codon_exn_num = str(start_codon_dict[gene_trans][start_codon_idx][6])

    
   
    # if the start codon is in the same exon as the breakpoin

    # calculating lens of intermediate exons
    length_intermed_exons = 0
    for i1 in range(start_codon_idx + 1, index_of_bp_exon):
        length_intermed_exons += int(start_codon_dict[gene_trans][i1][4])
    #print('length_intermed_exon', length_intermed_exons)
    # getting start & stop based on strand
    if strand == "+":
        start_codon_exon_start = int(start_codon_pos_list[2])
        start_codon_exon_end = int(start_codon_pos_list[3])
    else:
        start_codon_exon_start = int(start_codon_pos_list[3])
        start_codon_exon_end = int(start_codon_pos_list[2])

    '''
    Check if the breakpoint lies in the untranslated region (UTR).
    For the negative strand, this checks if the start codon is upstream of the breakpoint and the breakpoint is within the exon containing the start codon.
    For the positive strand, this checks if the start codon is downstream of the breakpoint and the breakpoint is within the exon containing the start codon.
    If so, return "In_Untranslated_Region" as the frame status.
    '''

    if (strand == "-" and  (start_codon < int(bkpt_coordinate) >= start_codon_exon_start) and start_codon_exn_num == "E"+ str(bp_exon)) or (strand == "+" and (start_codon > int(bkpt_coordinate) <= start_codon_exon_start) and start_codon_exn_num == "E"+ str(bp_exon)):

        # print(f'\n\n++++++++++++ TEST CASES: gene_transcript: {gene_trans} bkpt_coordinate: {bkpt_coordinate} strand: {strand} and start_codon: {start_codon} bkpt_exon_start: {start_codon_exon_start} and start_codon_exn_num: {start_codon_exn_num} bp_exon: {bp_exon}, start_codon_list: {start_codon_list} start_codon_idx: {start_codon_list.index(str(start_codon))}++++++++++++++\n\n')
        return result_dict, "In_Untranslated_Region"
    



      
    # Calculate the final length based on the conditions
    if start_codon == start_codon_exon_start:
        final_len = start_codon_exn_len + length_intermed_exons + bkpt_exn_len


    elif (strand == "+" and start_codon > int(start_codon_exon_start)) or \
            (strand == "-" and start_codon > int(start_codon_exon_end)):
        
        if start_codon_exn_num == "E" + str(bp_exon):
            if strand == "-":
                start_codon_to_exn_end = abs(bkpt_coordinate - int(start_codon_exon_start)) + 1
            elif strand == "+":
                # Working with Pytest
                start_codon_to_exn_end = abs(start_codon - int(start_codon_exon_end)) + 1 

                #current#
                # start_codon_to_exn_end = abs(bkpt_coordinate - int(start_codon)) + 1
                ##start_codon_to_exn_end = abs(start_codon - int(start_codon_exon_end)) + 1
                #start_codon_to_exn_end = abs(bkpt_coordinate - int(start_codon_exon_end)) + 1
            
                #start_codon_to_exn_end = abs(bkpt_coordinate - int(start_codon_exon_start)) + 1
            #print('test:',gene_trans,bp_exon,bkpt_coordinate,start_codon_to_exn_end)
            final_len = start_codon_to_exn_end
        else:
            start_codon_to_exn_end = abs(start_codon - int(start_codon_exon_end)) + 1
            final_len = start_codon_to_exn_end + length_intermed_exons + bkpt_exn_len
    else:
        final_len = bkpt_exn_len + length_intermed_exons
    modulo = (final_len) % 3
    frame = (final_len) % 3 + 1  # Adding 1 to convert to a 1-based frame index

    # Debug
#    logger.debug(f'gene_transcript:{gene_trans}:{exon_sent}\tstart_codon:{start_codon}\tstart_codon_exon_start:{start_codon_exon_start}\tstart_codon_exon_end:{start_codon_exon_end}\textra_bases:{extra_bases}\tfinal_length:{final_len}\tmodulo:{modulo}\tframe:{frame}')
    result_dict = {
        'gene_transcript': f'{gene_trans}:{exon_sent}',
        'coord': bkpt_coordinate,
        'condition': condition,
        'extraBases': extra_bases,
        'final_length': final_len,
        'length_intermediate': length_intermed_exons, 
        'modulo': modulo,
        'frame': frame
    }
    results_list.append(result_dict)
    #df = pd.DataFrame(results_list)
    #df.to_csv('output.csv', index=False)
    return result_dict, frame


def find_3p_start_codon_loc(gene_trans3p, strand_3p, exon_3p,
                            exon_loc_3p, chr_3p, coord_3p):
    """
    Function to compute the location of the breakpoint within the exon of the
    3 prime gene partner.
    The function adds, subtracts a certain number of base pairs based on the location
    of the breakpoint within the exon consisting of the breakpoint in order to
    determine the coding length of the exon contributing towards the frame calculation
    for the gene partner

    Attributes
    ----------
    gene_trans_3p, strand_3p, exon_loc_3p chr_3p : str
    coord_3p, exon_3p  : int

    Variables
    ---------
    exn_st = exon start (Note:- This variable will be exon end on the negative strand)
    exn_end = exon end (Note:- This variable will be exon start on the negative strand)
    exn_len = length of the exon consisting of the breakpoint
    exn_num = Exon number within the transcript.
    len_ff = frame length that will be calculated in this function.
    extra_bases = Depending on the location of the breakpoint within the exon, these
                  many # of bases need to be added or subtracted from exon length

    Methods
    -------
    compute_frame:
    The method will take the following attributes
    1. exon_3p: exon of the 3p partner
    2. exon_name : Current exon from the gene:transcript dictionary
    3. additional bases to be added or subtracted depending on the location of the break
       length of the exon
    4. the gene:transcript combination to be queried from the SC file
    5. idx : The index of the exon

     The method returns the frame of the 3p end of gene involved in the fusion
     based on where the breakpt is located within the exon of the fusion.
    """

    coord_3p, exon_3p = int(coord_3p), str(exon_3p)    
    extra_bases = 0
    condition = "" 

    for idx, exon_name in enumerate(start_codon_dict[gene_trans3p]):
        exn_st, exn_end, exn_len, exn_num = int(exon_name[2]), int(exon_name[3]), int(exon_name[4]), str(exon_name[6])
        if strand_3p == "-":
                st_extra = abs(coord_3p - exn_end)
                end_extra = abs(coord_3p - exn_st)
        else:
                st_extra = abs(coord_3p - exn_st)
                end_extra = abs(coord_3p - exn_end)

  # 1 : 3p bkpt is at the END of the exon on the "-" strand or within the BODY
        if (exon_3p == exn_num and strand_3p == "-") and (exon_loc_3p == "End"):
            # 1A: If 3p is at the END, exon len is excluded with 1 extra base
            if coord_3p == exn_st:
                condition = "1A-End"
                extra_bases = 1
            # 1B: If 3p bkpt is within the body upstream of exon END (greater than col2 in SC file)
            if coord_3p > exn_st:
                condition = "1B-End"
                extra_bases = end_extra + 1            
            # 1C: If 3p bkpt is downtream of the exon END (less than the col2 - Intron Retention type)
            if coord_3p < exn_st:
                condition = "1C-End"
                return "NA"
        
        # 2 : 3p bkpt located at the END of the exon on the "+" strand or within the BODY
        if (exon_3p == exn_num and strand_3p == "+") and (exon_loc_3p == "End"):
            # 2A: If 3p bkpt is located exactly at the END
            if coord_3p == exn_end:
                condition = "2A+End"
                extra_bases = 1
            # 2B: If 3p bkpt is located before the END of the exon
            if coord_3p < exn_end:
                condition = "2B+End"
                extra_bases = end_extra + 1                
            # 2C: If a 3p bkpt is located downstream of exon END (Intron retention type case)
            if coord_3p > exn_end:
                condition = "2C+End"
                return "NA"
        
        # 3: 3p bkpt located at the START of the exon on the "-" strand or within the BODY
        if (exon_3p == exn_num and strand_3p == "-") and (exon_loc_3p == "Start"):
            # 3A: If 3p bkpt is exactly at the START of the exon; exon len is excluded with 1 extra base
            if coord_3p == exn_end:
                condition = "3A-Start"
                extra_bases = 1                            
             # 3B: If 3p bkpt is greater than exon START
            if coord_3p > exn_end:
                condition = "3B-Start"
                extra_bases = st_extra + 1
            # 3C: If 3p bkpt is less than exon START
            if coord_3p < exn_end:
                condition = "3C-Start"
                extra_bases = end_extra + 1
        
        # 4:3p bkpt located at the START of the exon on the "+" strand, within the BODY
        if (exon_3p == exn_num and strand_3p == "+") and (exon_loc_3p == "Start"): 
            # 4A: If 3p bkpt is exactly at the START of the exon (count only 1 bp from the exon)
            if coord_3p == exn_st:
                condition = "4A+Start"
                extra_bases = 1
            # 4B: If 3p bkpt is before the START of the exon; Intron retention type case
            if coord_3p < exn_st:
                condition = "4B+Start"
                return "NA"
            # 4C: If the 3p bkpt is after the START of the exon
            if coord_3p > exn_st:
                condition = "4C+Start"
                extra_bases = st_extra + 1

        #BODY and Intron conditions
        # 1 : 3p bkpt is at the END of the exon on the "-" strand or within the BODY
        if (exon_3p == exn_num and strand_3p == "-") and (exon_loc_3p == "Body"):
            # 1A/5A: If 3p bkpt is before the exon
            if coord_3p > exn_st and coord_3p < exn_end:
                condition = "5A-BodyNeg"
                extra_bases = st_extra + 1
        #Body Pos
        if (exon_3p == exn_num and strand_3p == "+") and (exon_loc_3p == "Body"):
            # 2B/6A: If 3p bkpt is located before the END of the exon
            if coord_3p < exn_end and coord_3p > exn_st:
                condition = "6A+BodyPos"
                extra_bases = end_extra + 1

    # Check if any condition is met and update the 'extra_bases' accordingly
        if condition:
            break
    if condition:
        return compute_frame(exon_3p, exon_name, extra_bases, gene_trans3p, idx, coord_3p, condition)
    else:
        return None, None





def find_5p_start_codon_loc(gene_trans5p, strand_5p, exon_5p,
                            exon_loc_5p, chr_5p, coord_5p):
    """
    Function to compute the location of the breakpoint within the exon of the
    3 prime gene partner.
    The function adds, subtracts a certain number of base pairs based on the
    breakpoint locationwithin the exon consisting of the breakpoint in order to
    determine the coding length of the exon contributing towards the frame
    calculation for the gene partner

    Attributes
    ----------
    gene_trans_5p, strand_5p, exon_loc_5p chr_5p : str
    coord_5p, exon_5p  : int

    Variables
    ---------
    exn_st = exon start (Note:- This variable will be exon end on the '-' strand)
    exn_end = exon end (Note:- This variable will be exon start on the '-' strand)
    exn_len = length of the exon consisting of the breakpoint
    exn_num = Exon number within the transcript.
    len_ff = frame length that will be calculated in this function.
    extra_bases = Depending on the location of the breakpoint within the exon,
                  these extra bases need to be added or subtracted from exon length

    Methods
    -------
    compute_frame:
    The method will take the following attributes
    1. exon_5p: exon of the 5p partner
    2. exon_name : Current exon from the gene:transcript dictionary
    3. extra bases to be added or subtracted depending on the location of the break
       length of the exon
    4. the gene:transcript combination to be queried from the SC file
    5. idx : The index of the exon

     The method returns the frame of the 5p end of gene involved in the fusion
     based on where the breakpt is located within the exon of the fusion.
    #TODO
    K change it to exon name
    start_i , end_i change it to extra bases and exn length
    gene_trans5p change it gene_trans
    st_extra,chang
    """

    coord_5p, exon_5p = int(coord_5p), str(exon_5p)    
    extra_bases = 0
    condition = ""
    exn_len = 0
    # Initialize a flag variable
    condition_met = False
    for idx, exon_name in enumerate(start_codon_dict[gene_trans5p]):
        exn_st, exn_end, exn_len, exn_num = int(exon_name[2]), int(exon_name[3]), int(exon_name[4]), str(exon_name[6])
        if strand_5p == "-":
                st_extra = abs(coord_5p - exn_end)
                end_extra = abs(coord_5p - exn_st)
        else:
                st_extra = abs(coord_5p - exn_st)
                end_extra = abs(coord_5p - exn_end)
        
        # 1 : 5p bkpt located at the END of the exon and "-" strand or within BODY
        if (exon_5p == exn_num and strand_5p == "-") and (exon_loc_5p == "End"):
            # 1A: If 5p breakpoint is exactly at the END
            if coord_5p == exn_st:
                condition = "1A-End"
                extra_bases = exn_len
            # 1B: If the breakpoint is located before the end of the exon
            if coord_5p > exn_st:
                condition = "1B-End"
                extra_bases = st_extra + 1
            # 1C: If the breakpoint is after the END of an exon
            if coord_5p < exn_st:
                condition = "1C-End"
                extra_bases = end_extra +1
                
        # 2: bkpt located at the END of the exon and "+" strand or within  BODY
        if (exon_5p == exn_num and strand_5p == "+") and (exon_loc_5p == "End"):
            # 2A: If breakpoint is EXACTLY at the END
            if coord_5p == exn_end:
                condition = "2A+End"
                extra_bases = exn_len
            # 2B: If breakpoint is BEFORE the END
            #this will satisfy if the breakpoint is less than 15, if its more will fall in 6A
            if coord_5p < exn_end:
                condition = "2B+End"
                if abs(coord_5p - exn_end) == 1:
                    extra_bases = exn_len - 1 # Include start_i = 1 in extra_bases hence exn_len -2 is changed 
                else:
                    exn_len = exn_len - (end_extra + 1)
                    extra_bases = exn_len 
            # 2C: If breakpoint is AFTER the end of the exon
            if coord_5p > exn_end:
                condition = "2C+End"
                extra_bases = end_extra + 1 
        #  3: bkpt is located at the START of the exon and "-" strand or within BODY
        if (exon_5p == exn_num and strand_5p == "-") and (exon_loc_5p == "Start" ):
            # 3A: If breakpoint is EXACTLY at the START of the exon
            if coord_5p == exn_end:
                condition = "3A-Start"
                extra_bases = 1
            # 3B: If 5p breakpoint is upstream of the exon start ; Very rare
            if coord_5p > exn_end:
                condition = "3B-Start"
                return "NA"
            # 3C: If bkpt is less than actual exon start, it lies within the body
            if coord_5p < exn_end:
                condition = "3C-Start"            
                extra_bases = st_extra + 1
                    
        # 4: 5p bkpt located at the START of the exon and '+' strand or within BODY
        if (exon_5p == exn_num and strand_5p == "+") and (exon_loc_5p == "Start"):
            # 4A: Exactly at the start of the exon
            if coord_5p == exn_st:
                condition = "4A+Start"
                extra_bases = 1
            # 4B: Before the start of the exon on the + strand ; very rare case
            if coord_5p < exn_st:
                condition = "4B+Start"
                return "NA"
            # 4C : bkpt is located AFTER the start of the exon on the '+' strand
            if coord_5p > exn_st:
                condition = "4C+Start"
                extra_bases = st_extra + 1

    #### BODY and Intron conditions
        if (exon_5p == exn_num and strand_5p == "-") and (exon_loc_5p == "Body"):           
            # 1B/5A: If the breakpoint is located before the end of the exon
            if coord_5p > exn_st and coord_5p < exn_end:
                condition = "5A-BodyNeg"
                extra_bases = st_extra + 1
        
        if (exon_5p == exn_num and strand_5p == "+") and (exon_loc_5p == "Body"):            
            # 2B/6A: If breakpoint is BEFORE the END
            if coord_5p < exn_end and coord_5p > exn_st:
                condition = "6A+BodyPos" # ******************* diff than 2B
                extra_bases = st_extra + 1

    # Check if any condition is met and update the 'extra_bases' accordingly
        if condition:
            break
    if condition:
        return compute_frame(exon_5p, exon_name, extra_bases, gene_trans5p, idx, coord_5p, condition)
    else:
       
        logger.debug(f'No condition met for gene {gene_trans5p} exon_5p, exon_name, extra_bases, gene_trans5p, idx, coord_5p, exon_5p: {exon_5p}, exon_name: {exon_name}, extra_bases: {extra_bases}, gene_trans5p: {gene_trans5p}, idx: {idx}, coord_5p: {coord_5p}, {condition}')
        return None, None

def parse_coordinates(header_index, coordinate, j):
    """
    The function parses the field consisting of the exon annotation
    to extract the gene , transcript, strand, exon , exon location
    and chromosomal coordinates of a given fusion breakpoint

    Attributes
    ----------
    header_index - index of column consisting of the exon annotation
    coordinate - chromsomal coordinate of the breakpoint
    j = The current row of the incoming fusion dataframe being iterated.

    """
    gene = j[header_index].split("|")[2]
    transcript = j[header_index].split("|")[3]
    gene_trans = str(gene) + ":" + str(transcript)
    strand = j[header_index].split("|")[0]
    exon = j[header_index].split("|")[1].split("_")[1]
    exon_loc = j[header_index].split("|")[1].split("_")[0]
    chrs = j[coordinate].split(":")[0]
    numerical_coord = int(j[coordinate].split(":")[1])
    return [gene, transcript, gene_trans, strand, exon,
            exon_loc, chrs, numerical_coord]


def main():
    parser = ArgumentParser()
    try:
        args = parser.parse_args()
    except:
        sys.exc_info()[1]
        logger.error(f"Please provide above given arguments")
        exit()

    # Read input aruguments
    input_file = args.i
    codon_file = args.c


def frame_calculation(input_fusions):
    """
    Function to calculate the frame of a given fusion candidates

    Attributes
    -----------
    input_fusions :- Data frame of incoming fusion candidates

    """

    # Getting resource path
    logger.info("-- estimating frame of fusion candidates --")

    debug_df = pd.DataFrame()
    # Reading the codon file as well as the input file
    fusion_file = utils.read_df_or_file(input_fusions)
    start_codon_file = utils.read_df_or_file(config["start_codon_file"])

    global start_codon_dict
    start_codon_dict = {}
    for i in start_codon_file:
        i = i.strip().split("\t")
        gene_transcript_combo = str(i[7]) + ":" + str(i[8])
        curr_value = i[0:7]
        if start_codon_dict.get(gene_transcript_combo) is not None:
            start_codon_dict[gene_transcript_combo].append(curr_value)
        else:
            start_codon_dict[gene_transcript_combo] = [curr_value]

    fus_header = ["Frame_5p", "Frame_3p", "Frame_Status"]
    final_list = []
    fus2 = 0
    #bar = progressbar.ProgressBar(max_value=len(input_fusions))
    for idx, line in enumerate(tqdm(fusion_file,  leave=False, ascii="█▖▘▝▗▚▞=")):
        frame_5p, frame_3p = "NA", "NA"
        j = line.strip().split("\t")
        if idx == 0:
            fus_header = j + fus_header
            header_5p_index = j.index("5'_Exon_Annotation")
            coordinate_5p = j.index("5'co-ordinate")
            header_3p_index = j.index("3'_Exon_Annotation")
            coordinate_3p = j.index("3'co-ordinate")

        else:
            fus = j[0]
            if j[header_5p_index] != "NA":
                info_list_5p = parse_coordinates(header_5p_index, coordinate_5p, j)
                #(len(find_5p_start_codon_loc(*info_list_5p[2:])), find_5p_start_codon_loc(*info_list_5p[2:]))
                results_5p, frame_5p = find_5p_start_codon_loc(*info_list_5p[2:])
            if j[header_3p_index] != "NA":
                info_list_3p = parse_coordinates(header_3p_index, coordinate_3p, j)
                results_3p, frame_3p = find_3p_start_codon_loc(*info_list_3p[2:])
            frame_status = get_frame_status(frame_5p, frame_3p)
            
            if logger.isEnabledFor(logging.DEBUG):
                df1 = pd.DataFrame([results_5p]).add_suffix('_5p')
                df2 = pd.DataFrame([results_3p]).add_suffix('_3p')
                combined = pd.concat([df1, df2], axis=1)
                debug_df = pd.concat([debug_df, combined])  
                debug_df['frame_status'] = frame_status


            temp_row = [frame_5p, frame_3p, frame_status]
            new_j = j + temp_row
            final_list.append(new_j)
            fus2 = fus
        #bar.update(idx)

    final_df = pd.DataFrame(final_list, columns=fus_header)

    if logger.isEnabledFor(logging.DEBUG):
        debug_df = debug_df.reset_index()
        debug_out_file = os.path.join(utils.out_path , "DEBUG_frame_calculation.tsv")
        logger.debug(f"Debugging info for frame calculation saved to {debug_out_file}")
        debug_df.to_csv(debug_out_file, sep="\t", index=False)
    return final_df


if __name__ == "__main__":
    main()
