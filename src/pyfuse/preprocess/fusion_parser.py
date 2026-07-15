#!/usr/bin/env python3
"""A parser for various fusion caller output formats to extract breakpoints"""
__status__ = "Prototype"
__date__ = "01/26/2023"

import sys
import logging
import numpy as np
import pandas as pd
import pybedtools as pbt
from pyfuse.utils.common_utils import utils, config
pd.options.mode.chained_assignment = None

logger = logging.getLogger(__name__)


class ParseFusion():
    """
    A class to preprocess input fusion breakpoints

    Attributes
    ----------
    fusion_bkpt_content: DF or file stream
        path of input breakpoint file
    res_path: path
        path of resource files needed for package annotation
    format: str
        fusion caller name to preprocess input breakpoint file according to
        respective caller output format

    Methods
    -------
    preprocess_bkpt(object):
        handles all the preprocessing
    filter_bkpts_with_bed(object):
        checks if input brakpoints are within the input bed region return boolean
    """

    def __init__(self, preprocess_instance):
        self.preprocess = preprocess_instance

    def parse_fusion(self):
        method_name = f"parse_{ self.preprocess.format}_bkpt"

        logger.info(f"User provided input format/caller is: {self.preprocess.format}")

        # Check if the method exists
        if hasattr(self, method_name):
            # Dynamically get the method and call it
            logger.info(f"Parsing { self.preprocess.format} output")
            method = getattr(self, method_name)
            method()
        else:
            logger.error(f"Unsupported format: {self.preprocess.format}")

    def parse_star_bkpt(self):
        # parsing 5' and 3' breakpoints from Star fusion output
        df = pd.read_csv(self.preprocess.fusion_bkpt, sep='\t')
        try:
            cols5, cols3 = ['chr5', 'coor5', 'strand5'], ['chr3', 'coor3', 'strand3']
            bkpt_df = df[['LeftBreakpoint', 'RightBreakpoint']]
            bkpt_df[cols5] = bkpt_df['LeftBreakpoint'].str.split(":", expand=True)
            bkpt_df[cols3] = bkpt_df['RightBreakpoint'].str.split(":", expand=True)
            self.preprocess.main_bkpt_df = bkpt_df
        except KeyError as e:
            sys.exit(f'Error: {e}\n please check if the caller format used in --input_format is correct or check if the input file has the required columns for Star fusion parsing')

    def parse_default_bkpt(self):

        try:
            df = pd.read_csv(self.preprocess.fusion_bkpt, sep='\t')
            self.preprocess.main_bkpt_df = df.iloc[:, :4]
            self.preprocess.main_bkpt_df.columns = ['chr5', 'coor5', 'chr3', 'coor3']

            #if first column first row doesn't start with "chr" then raise error
            if not self.preprocess.main_bkpt_df.iloc[0, 0].startswith("chr"):
                raise ValueError("The first column of the input file does not start with 'chr'. Please check the input file format or use the correct --input_format option.")
        except KeyError as e:
            sys.exit(f'Error: {e}\n please check if the caller format used in --input_format is correct or check if the input file has the required columns for default parsing')
        except ValueError as e:
            sys.exit(f'Error: {e}\n please check if the caller format used in --input_format is correct or check if the input file has the four required columns for default parsing')

    def parse_arriba_bkpt(self):
        try:
            df = pd.read_csv(self.preprocess.fusion_bkpt, sep='\t')
            df = df[['breakpoint1', 'breakpoint2']]
            df[['chr5', 'coor5']] = df['breakpoint1'].str.split(':', 1, expand=True)
            df[['chr3', 'coor3']] = df['breakpoint2'].str.split(':', 1, expand=True)
            self.preprocess.main_bkpt_df = df
        except KeyError as e:
            sys.exit(f'Error: {e}\n please check if the caller format used in --input_format is correct or check if the input file has the required columns for Arriba parsing')

    def parse_tophat_bkpt(self):
        try:
            df = pd.read_csv(self.preprocess.fusion_bkpt, sep='\t', header=None)
            self.preprocess.main_bkpt_df = df.iloc[:, [2, 3, 5, 6]]
            self.preprocess.main_bkpt_df.columns = ['chr5', 'coor5', 'chr3', 'coor3']
        except KeyError as e:
            sys.exit(f'Error: {e}\n please check if the caller format used in --input_format is correct or check if the input file has the required columns for TopHat parsing')

    def parse_fusion_catcher_bkpt(self):
        try:
            df = pd.read_csv(self.preprocess.fusion_bkpt, sep='\t')
            df[['breakpoint1', 'breakpoint2']] = df[['Fusion_point_for_gene_1(5end_fusion_partner)', 'Fusion_point_for_gene_2(3end_fusion_partner)']]
            df[['chr5', 'coor5']] = df['breakpoint1'].str.split(':|:', expand=True)[[0, 1]]
            df[['chr3', 'coor3']] = df['breakpoint2'].str.split(':|:', expand=True)[[0, 1]]
            self.preprocess.main_bkpt_df = df[['chr5', 'coor5', 'chr3', 'coor3']]
        except KeyError as e:
            sys.exit(f'Error: {e}\n please check if the caller format used in --input_format is correct or check if the input file has the required columns for FusionCatcher parsing')

    def parse_longgf_bkpt(self):
        try:
            df = pd.read_csv(self.preprocess.fusion_bkpt, sep='\t', header=None)
            df[['breakpoint1', 'breakpoint2']] = df.iloc[:, 1].str.split(' ', expand=True).iloc[:, [7, 9]]
            df[['chr5', 'coor5']] = df['breakpoint1'].str.split(':', 2, expand=True)
            df[['chr3', 'coor3']] = df['breakpoint2'].str.split(':', 2, expand=True)
            self.preprocess.main_bkpt_df = df[['chr5', 'coor5', 'chr3', 'coor3']]
        except KeyError as e:
            sys.exit(f'Error: {e}\n please check if the caller format used in --input_format is correct or check if the input file has the required columns for LongGF parsing')

    def parse_fusion_inspector_bkpt(self):
        try:
            df = pd.read_csv(self.preprocess.fusion_bkpt, sep='\t')
            df[['chr5', 'coor5']] = df['LeftBreakpoint'].str.split(':|:', expand=True)[[0, 1]]
            df[['chr3', 'coor3']] = df['RightBreakpoint'].str.split(':|:', expand=True)[[0, 1]]
            self.preprocess.main_bkpt_df = df[['chr5', 'coor5', 'chr3', 'coor3']]
        except KeyError as e:
            sys.exit(f'Error: {e}\n please check if the caller format used in --input_format is correct or check if the input file has the required columns for Fusion Inspector parsing')
