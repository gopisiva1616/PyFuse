#!/usr/bin/env python3


'''Overlap fusion candidates from fusion detection pipelines with fusions from COSMIC and GTEX'''



import os
import sys
import csv
import time
import logging
import argparse
import pandas as pd
from pyfuse.utils.common_utils import utils, config

start_time = time.time()
logger = logging.getLogger(__name__)
csv.field_size_limit(sys.maxsize)


def validatefilepaths(path):
    if os.path.exists(path) is False:
        raise Exception("ERROR: The given path" + path + "does not exist \
                        or the path is incorrect \n Please check the path")
    if os.path.getsize(path) == 0:
        raise Exception("ERROR: The given path" + path + "does not exist or is inaccessible")
    return path


def argumentparser():
    help_text = '''The script annotates incoming fusion candidates if they are
                present with GTEX and COSMIC databases and assigns additional
                coloumns to the input file accordingly'''
    parser = argparse.ArgumentParser(description=help_text, epilog=help_text)
    parser.add_argument("-i", help="Path to input tab delimited fusion file with first column as fusion name \
                        example, EML4-ALK", required="True", type=validatefilepaths)
    parser.add_argument("-c", help="COSMIC fusion file")
    parser.add_argument("-g", help="GTEX fusion file")
    parser.add_argument("-o", help="Path to output directory", required="True",
                        type=validatefilepaths)
    return parser


def main():

    parser = argumentparser()
    try:
        args = parser.parse_args()
    except:
        sys.exc_info()[1]
        print("Please provide above given arguments")
        exit()

    # Read input aruguments
    input_file = args.i
    cosmic_file = args.c
    gtex_file = args.g


def annotate_gtex_cosmic(input_fusions):
    """
    Function annotate incoming fusion candidate if it present within
    the database of GTEX and COSMIC fusions

    Attribute
    ---------
    input_fusins: dataframe of incoming fusions
    """

    # getting resource path
    logger.info("-- Annotating fusions with GTEx and COSMIC data")
    res_path = utils.get_res_path()

    # Reading the gtex and cosmic files
    gtex_fusions = utils.read_df_or_file(config['gtex_fusions'])
    cosmic_fusions = utils.read_df_or_file(config['cosmic_fusions'])
    fusions_file = utils.read_df_or_file(input_fusions)

    # Make a dictionary of cosmic fusions
    cosmic_fusions_dict = {}
    for idx, line in enumerate(cosmic_fusions):
        if idx != 0:
            i = line.strip().split('\t')
            if i[0] not in cosmic_fusions_dict.keys():
                cosmic_fusions_dict[str(i[0])] = str(i[1])

    # Make a dictionary of gtex fusions
    gtex_fusions_dict ={}
    for idx, line in enumerate(gtex_fusions):
        if idx != 0:
            i = line.strip().split('\t')
            if i[0] not in gtex_fusions_dict.keys():
                gtex_fusions_dict[str(i[0])] = [i[1], i[2], i[3]]

    # Make a header for the resulting output file
    fus_header = ["Present_in_COSMIC", "Histology", "Present_in_GTEX",
                  "Average_Expression", "Number_of_Tissues_that_contain_fusion",
                  "Tissue_Names"]

    final_list = []
    for idx, line in enumerate(fusions_file):
        i = line.strip().split('\t')
        if idx == 0:
            fus_header = i + fus_header
            continue

        five_p = i[1].split("-")[0]
        three_p = i[1].split("-")[1]
        fus_can = str(five_p)+"-"+str(three_p)
        if fus_can in cosmic_fusions_dict.keys():
            cosmic_bin = "Yes"
            temp_row = [cosmic_bin, cosmic_fusions_dict[fus_can]]
        else:
            cosmic_bin = "No"
            temp_row = [cosmic_bin, "NA"]

        if fus_can in gtex_fusions_dict.keys():
            gtex_bin = "Yes"
            new_temp = [gtex_bin, gtex_fusions_dict[fus_can][0],
                        gtex_fusions_dict[fus_can][1],
                        gtex_fusions_dict[fus_can][2]]
            temp_row.extend(new_temp)
        else:
            gtex_bin = "No"
            new_temp = [gtex_bin, "NA", "NA", "NA"]
            temp_row.extend(new_temp)

        new_i = i + temp_row
        final_list.append(new_i)

    final_df = pd.DataFrame(final_list, columns=fus_header)
    #final_df.to_csv('gtext_inte.txt', header=True)
    return final_df


if __name__ == '__main__':
    main()
