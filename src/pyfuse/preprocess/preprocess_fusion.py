#!/usr/bin/env python3

"""A class to preprocess input fusion breakpoints"""
__status__ = "Prototype"
__date__ = "05/05/2020"

import os
import sys
import time
import logging
import progressbar
import numpy as np
import pandas as pd
import pybedtools as pbt
from .fusion_parser import ParseFusion
from pyfuse.utils.common_utils import utils, config
pd.options.mode.chained_assignment = None

logger = logging.getLogger(__name__)


class PrepareInputBkpt():
    """
    A class to preprocess input fusion breakpoints

    Attributes
    ----------
    fusion_bkpt_content: DF or file stream
        path of input breakpoint file
    out_path: path
        path to write output and log files
    res_path: path
        path of resource files needed for package annotation
    format: str
        fusion caller name to preprocess input breakpoint file according to
        respective caller output format
    target_df: dataframe
        bed DF to filter breakpoint that are within this bed region
    input_fusion_count: int
        count of input breakpoints
    exon_filt_count: int
        count of breakpoint after filtering with exon-region bed file
    target_filt_count: int
        count of breakpoint after filtering with target region bed file
    excluded_bkpt: list
        breakpoint that were filtered with exon-regions or target bed files

    Methods
    -------
    preprocess_bkpt(object):
        handles all the preprocessing
    filter_bkpts_with_bed(object):
        checks if input brakpoints are within the input bed region return boolean
    """

    def __init__(self, fusion_input_file, format, target_df, out_path):

        self.target_df = target_df
        self.final_bkpt_list = []
        self.fusion_bkpt = fusion_input_file
        self.format = format
        self.excluded_bkpt = pd.DataFrame()
        self.all_gene_df = utils.read_df_or_file(config['all_gene_exons'])

    def preprocess_bkpt(self):
        '''
        Parameters
        ----------
        input_fusion_count: int
            count of input brakpoints
        exon_filt_count: int
            count of breakpoint after filtering with exon-region bed file
        target_filt_count: int
            count of breakpoint after filtering with target region bed file
        excluded_bkpt: list
            breakpoint that were filtered with exon-regions or target bed files
        final_bkpt_list: list
            list of final breakpoints after preprocessing

        methods and functions
        --------------------
        utils.extract_fusion_bkpt(str):
            external function to extarct breakpoints from a input breakpoint line
        filter_bkpts_with_bed():
            checks if input brakpoints are within the input bed region return boolean
        '''
        start = time.process_time()
        input_fusion_count, exon_filt_count = 0, 0

        # parse input bkpts acording to input file format
        parser = ParseFusion(self)
        # Call the parse_file method, which will update self.bkpt
        parser.parse_fusion()

        # convert input bkpts to bedpe format for filtering
        self.main_bkpt_df['Fusion_id'] = np.arange(len(self.main_bkpt_df)) + 1
        input_fusion_count = len(self.main_bkpt_df)
        logger.info(f'Total input fusions: {input_fusion_count}')
        self.prep_df_for_intersect(self.main_bkpt_df.copy())

        # ----- filter bkpt further if target bed is provided ----
        if self.target_df.empty:
            logger.info(f'Target filtered: target bed not provided')
        else:
            logger.info(f'Target bed provided as input')
            target_bed = utils.df2pbt(self.target_df)
            intersect_df = self.intersect_bedpe(self.bkpt2filt_bed, target_bed, 'either')
            overlap, nonoverlap = self.extract_overlap_bkpts(self.main_bkpt_df,
                                                             intersect_df,
                                                             'non_target_bkpt')
            self.excluded_bkpt = nonoverlap
            self.main_bkpt_df = overlap
            self.prep_df_for_intersect(self.main_bkpt_df)
            logger.info(f'Non target fusions: {len(nonoverlap)}; Target fusions: '
                        f'{len(overlap)}')

        # ------- Gene exon filter -------
        # Filter bkpts using all gene coordinates to ensure bkpts of withins genes

        all_genes_bed = utils.df2pbt(utils.slop_bed(self.all_gene_df, config['gene_filter_padding']))
        intersect_df2 = self.intersect_bedpe(self.bkpt2filt_bed, all_genes_bed, 'both')

        if intersect_df2.empty:
            logger.error('None of the input breakpoints are part of the inbuilt gene models for adding gene annotation. You can try with custom reference files')
            sys.exit()
        overlap2, nonoverlap2 = self.extract_overlap_bkpts(self.main_bkpt_df,
                                                           intersect_df2,
                                                           'outside genes')
        self.main_bkpt_df = overlap2
        exon_filt_count = len(nonoverlap2)
        self.excluded_bkpt = nonoverlap2 if self.excluded_bkpt.empty else \
            self.excluded_bkpt.append(nonoverlap2)

        # resetting fusion-ids from 1
        self.main_bkpt_df['Fusion_id'] = np.arange(len(self.main_bkpt_df)) + 1
        logger.info(f'fusions outside genes: {exon_filt_count}; Within gene: '
                    f'{len(overlap2)}')
        logger.info(f'Preprocessed input fusion breakpoints')
        logger.info(f'Breakpoints to be annotated: {len(self.main_bkpt_df)}')
        logger.info(f'Pre-pocessing time: {(time.process_time() - start)/60 :.2} mins')
        return self.main_bkpt_df, self.excluded_bkpt, input_fusion_count

    def prep_df_for_intersect(self, df):
        self.bkpt2filt_df = utils.bed2bedpe(df, 1, 1)
        self.bkpt2filt_df.drop(['coor5', 'coor3'], axis=1, inplace=True)
        self.bkpt2filt_bed = utils.df2pbt(self.bkpt2filt_df)

    @staticmethod
    def intersect_bedpe(bkpt_df, main_df, type_):
        try:
            #main_df.to_dataframe().to_csv('all_genes.tsv', sep='\t')
            return bkpt_df.pair_to_bed(main_df, type=type_).to_dataframe()
        except pd.errors.EmptyDataError:
            logger.error('error occured during bkpt to bed conversion ')
            return pd.DataFrame(columns=bkpt_df.to_dataframe().columns)

    @staticmethod
    def extract_overlap_bkpts(orig_df, intersect_df, filter_msg):
        cols = ['chr5', 'coor5', 'chr3', 'coor3', 'Fusion_id']
        overlaps = orig_df[orig_df["Fusion_id"].isin(set(intersect_df['thickStart']))][cols]
        nonoverlaps = orig_df[~orig_df["Fusion_id"].isin(set(intersect_df['thickStart']))][cols]
        nonoverlaps['FILTER_REASON'] = filter_msg
        return overlaps, nonoverlaps
