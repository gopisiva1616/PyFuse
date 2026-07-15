#!/usr/bin/env python3

"""Common utility functions for PyFuse project"""
__status__ = "Prototype"
__date__ = "05/05/2020"

from importlib.resources import as_file, files
import io
import os
import sys
import yaml
import logging
import re
import numpy as np
import pandas as pd
import pybedtools as pbt
from pathlib import Path

logger = logging.getLogger(__name__)

class Utils():

    def __init__(self):
        #self.res_path  = self.get_res_path()
        self.config = self.get_config_vars()
        

    def read_df_or_file(self, input, path=None, df=False):
        """
        identify if input is dataframe or a file
        and return a sting buffer that is loopable
        line by line
        """
        def read_df(input):
            rstream = io.StringIO()
            input.to_csv(rstream, sep='\t', index=False)
            rstream.seek(0)
            return rstream
    

        if isinstance(input, pd.DataFrame):
            rstream = read_df(input)
        elif os.path.isfile(input) and input.endswith('.gz'):
            rstream = read_df(self.read_bed(input))
        elif os.path.isfile(input):
            logger.info(f'Reading file : {input}')
            with open(input, 'r') as f:
                rstream = f.readlines()
        else:
            logger.error(f'Invalid file: {input}')
        
        if df:
            rstream.seek(0)
            return pd.read_csv(rstream, sep='\t')
        return rstream

    @staticmethod
    def gene_pair_to_html_links(pair):

        try:
            gene1, gene2 = pair.split('-')
            config = utils.get_config_vars()
            url1 = f'{config["genecards_url"]}{gene1}'
            url2 = f'{config["genecards_url"]}{gene2}'
            return f'<a href="{url1}" target="_blank">{gene1}</a>-<a href="{url2}" target="_blank">{gene2}</a>'
        except ValueError as e:
            logger.warning(f'Unable to add  genecard link for the gene-pair: {pair} in the html report')
            return pair

    @staticmethod
    def coord_to_ucsc_link(coord):
        try:
            chrom, pos = str(coord).rsplit(':', 1)
            pos_int = int(pos)
            config = utils.get_config_vars()
            url = f'{config["ucsc_url"]}{chrom}:{pos_int - 50}-{pos_int + 50}'
            return f'<a href="{url}" target="_blank">{coord}</a>'
        except (ValueError, KeyError):
            logger.warning(f'Unable to add UCSC link for coordinate: {coord} in the html report')
            return coord

    @staticmethod
    def exon_annotation_to_refseq_link(annotation):
        """Add NCBI link to RefSeq transcript in 'strand|loc|gene|transcript' strings."""
        if not isinstance(annotation, str):
            return annotation

        parts = annotation.split('|')
        if len(parts) != 4:
            return annotation

        transcript = parts[-1].strip()
        if transcript in {'.', ''}:
            return annotation

        if not re.match(r'^[NX][MR]_\d+(?:\.\d+)?$', transcript):
            return annotation

        url = f'https://www.ncbi.nlm.nih.gov/nuccore/{transcript}'
        parts[-1] = f'<a href="{url}" target="_blank">{transcript}</a>'
        return '|'.join(parts)
            


    def check_path(self, path=None, msg=None):
        """
        checks the existence of input file and raise
        error with custom or default message
        """
        if not os.path.exists(path):
            if not msg:
                logger.exception(f'Required file/folder does not exist in {path}')
            else:
                logger.exception(msg)
            sys.exit()
        else:
            return os.path.abspath(path)



    def verify_res_files(self, res_path):
        """
        verifies if the input resource folder has all the required
        files for the function of the tool.
        """
        missing_files = []
        logger.info('checking required resource files')
        req_resource_files = self.get_config_vars()['required_res_files']
        res_files = os.listdir(res_path)

        for f in req_resource_files:
            if f not in res_files:
                missing_files.append(f)

        if missing_files:
            sys.exit(logger.error(f"Required resource files {missing_files} "
                                  f" does not exist in {res_path}"))

    def read_bed(self, file):
        ''' reading and returning exon_region resource file '''
        return pd.read_csv(file, sep='\t')

    @staticmethod
    def load_yaml(path: Path) -> dict:
        return yaml.safe_load(path.read_text(encoding="utf-8")) or {}

    def get_config_vars(self, user_config_path: str | None = None) -> dict:
        defaults_txt = (files("pyfuse.config") / "settings.yaml").read_text(encoding="utf-8")
        cfg = yaml.safe_load(defaults_txt) or {}

        if user_config_path:
            user_cfg = self.load_yaml(Path(user_config_path).expanduser())
            cfg.update(user_cfg)  # shallow merge; can do deep merge if needed
        return cfg

    @staticmethod
    def bed2bedpe(df, p1, p2):
        df[['coor5', 'coor3']] = df[['coor5', 'coor3']].astype(int)
        df['start5'], df['end5'] = df['coor5'] - p1, df['coor5'] + p2
        df['start3'], df['end3'] = df['coor3'] - p1, df['coor3'] + p2
        df = df[['chr5', 'start5', 'end5', 'chr3', 'start3', 'end3', 'coor5',
                 'coor3', 'Fusion_id']]
        return df

    def intersect_bed_df(self, bkpt_df, main_df):
        try:
            bkpt_bed = self.df2pbt(bkpt_df)
            main_bed = self.df2pbt(main_df)
            cols = list(main_df.columns) + list(bkpt_df.columns)
            b = (main_bed.intersect(bkpt_bed, wb=True).to_dataframe(header=None,
                 disable_auto_names=True))
            b.columns = cols
            b.columns = b.columns.str.replace(r'\d+|_x|_y', '', regex=True)
            b = b.loc[:, ~b.columns.duplicated()]
            return b
        except:
            logger.warning(sys.exc_info())
            return pd.DataFrame()

    @staticmethod
    def df2pbt(df):
        if isinstance(df, pd.DataFrame):
            return pbt.BedTool.from_dataframe(df)
        return df

    @staticmethod
    def slop_bed(bed_df, padding):
        if isinstance(bed_df, io.StringIO):
            bed_df.seek(0)
            bed_df = pd.read_csv(bed_df, sep='\t')

        bed_df['start'] = np.maximum(bed_df['start'].astype(int) - padding, 1)
        bed_df['end'] = bed_df['end'].astype(int) + padding
        return bed_df

    @staticmethod
    def extract_fusion_bkpt(line):
        """ Extracts chrom and co-ordinates of two breakpoint from string"""
        line_split = line.strip().split('\t')
        if line_split[0].startswith('chr'):
            bkpt1 = [line_split[0], int(line_split[1])]
            bkpt2 = [line_split[2], int(line_split[3])]
            return bkpt1, bkpt2


utils = Utils()
config = utils.config