#!/usr/bin/env python3

"""Annotate fusions as artifacts based on a database of black list genes"""

__status__ = "Development"

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
    help_text = "Annotate a given fusion for black list genes"
    parser = argparse.argumentparser(description=help_text, epilog=help_text)
    parser.add_argument("-i", help="input fusion file with first column as fusion name \
                        example, EML4-ALK", required="True", type=ValidatefilePaths)
    parser.add_argument("-b", help="List of black list fusions", type=ValidatefilePaths)
    return parser


def main():
    parser = argumentparser()
    try:
        args = parser.parse_args()
    except:
        sys.exc_info()[1]
        print(f'Please provide above given arguments')
        exit()

    input_file = args.i
    bl_list_file = args.b

def annotate_bl_list(input_fusions):
    """
    Annotate gene fusion candidates for their presence 
    within black list genes
    """
    logger.info("-- Adding Annotation for Black list")
    fusions_file = utils.read_df_or_file(input_fusions)
    black_list = utils.read_df_or_file(config['black_list'])
    header = ["General_Information", "\t", "\t", "\t", "\t",\
            "\t","Exon_Information", "\t", "\t", "\t", "\t",\
            "Read_Information", "\t", "\t", "Additional_Information"]
    bl = []
    fusions = []
    
    for i,item in enumerate(black_list):
        if i == 0:
            header_line = item.strip().split('\t')
        else:
            item = item.strip().split('\t')
            item_new = [item[0].replace("~", "-"),item[1]]
            bl.append(item_new)
    
    bl_df = pd.DataFrame(bl, columns = header_line)
    
    for i,item in enumerate(fusions_file):
        if i == 0:
            header = item.strip().split("\t")
        else:
            fusions.append(item.strip().split("\t"))

    fusions_df = pd.DataFrame(fusions , columns = header)
    
    #mergedDf = pd.merge(fusions_df, bl_df)
    mergedDf = pd.merge(fusions_df, bl_df, how='left',on="5'-3'Gene_Partners")
    mergedDf.fillna('NA', inplace=True)
    #final_df.to_csv('gtext_inte.txt', header=True)


    return mergedDf 



if __name__ == '__main__':
    main()
