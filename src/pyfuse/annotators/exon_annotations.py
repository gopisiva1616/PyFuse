#!/usr/bin/env python3

"""Annotate fusion breakpoints with exon-level information"""
__status__ = "development"
__date__ = "05/05/2020"

import os
import re
import time
import logging
import warnings
import argparse
import itertools
import numpy as np
import pandas as pd
from tqdm import tqdm
import pybedtools as pbt
from pyfuse.utils.common_utils import utils, config

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-a', dest='fusion_bkpt', required=True,
                        help="tsv file with SVs from happy", type=os.path.abspath)
    parser.add_argument('-n', dest='out_name', help="output name")
    parser.add_argument('-o', dest='out_path', type=os.path.abspath,
                        help="intersected output vcf")
    args = parser.parse_args()

    final_annoated_df = annotate_fusion_exons(args.fusion_bkpt, utils.res_path)
    if args.out_path:
        if args.out_name:
            out_file = args.out_path + '/' + args.out_name + '.tsv'
        else:
            out_file = args.out_path + '/fusion_exon_annoatations.tsv'
        #final_annoated_df.to_csv(out_file, sep='\t', index=False)


def annotate_fusion_exons(fusion_bkpt_df):
    '''
    Entry point for this script to read each input breakpoint and process
    them using dowstream functions

    Parameters
    ----------
    fusion_bkpt : list of bkpt1 and bkpt2 eg [chr:765342, chr:765342]
        5' and 3' fusion breakpoint location eg:chr:765342
    exon_region_df_dict : dict
        read from pkg config file, contains file name for start, stop and end, files
    main_df: Dataframe
        stores final output of each processed breakpoints

    Functions and classes:
    ---------------------
    df_from_csv(): function to read list of resource csv files that are in
    config's exon_region_dict and returns a dict of each file content as DF in a dict

    ProcessBreakpoints(): Main class that completely process each brakpint and
                          returns a df of annotation for input bkpts

    '''
    logger.info("-- Running exon-level annotation")
    main_list = []
    start = time.process_time()
    total_fusions = fusion_bkpt_df['Fusion_id'].max()
    bkpt5_df, bkpt3_df = prep_intersect_both_bkpts(fusion_bkpt_df)
    # bkpt5_df.to_csv('exon_annots5og.tsv', sep='\t', index=False)
    # bkpt3_df.to_csv('exon_annots3og.tsv', sep='\t', index=False)
    all_genes_df = utils.read_df_or_file(config['all_gene_exons'])
    annotator_obj = ProcessBreakpoints(bkpt5_df, bkpt3_df, all_genes_df)
    #bar = progressbar.ProgressBar(max_value=total_fusions)


    for fusion_id in tqdm(range(total_fusions), leave=False,
            ascii="█▖▘▝▗▚▞="):
        annotator_obj.run_annotation(fusion_id + 1)
        main_list.append(annotator_obj.fusion_df)
        #bar.update(fusion_id)
        
    main_df = pd.concat(main_list)
    # repalce intron annots like -|-|GPCPD1|- with NA
    main_df.replace(re.compile(r'-\|-\|.*\|-'), 'NA', inplace=True)
    #main_df['IGV link'] = main_df.apply(lambda row: get_igv_link(row), axis=1)
    logger.info(f'total exon annotation time: {round((time.process_time() - start)/60, 2)}mins')
    #main_df.to_csv('exon_annots.tsv', sep='\t', index=False)

    return main_df


def get_group(groupby, id):
    try:
        return groupby.get_group(id)
    except KeyError:
        return None


def prep_intersect_both_bkpts(fusion_bkpt_df):

    def process_df(bkpt_df, exons_df, all_genes_df):
        df = utils.intersect_bed_df(bkpt_df, exons_df)
        fdf = df
        df_na = bkpt_df[~bkpt_df.Fusion_id.isin(set(df.Fusion_id))]

        if not df_na.empty:
            ''' if not present in either exon start/end/body(exons_df), check if
            it is present in nearby intron(allgene_df (exon_df with 300bp padding)) '''

            df_na = utils.intersect_bed_df(df_na, all_genes_df)
            df_na.columns = df_na.columns.str.replace(r'\d+|_x|_y', '')

            df_na = df_na[['chr', 'coor', 'gene', 'Fusion_id']].drop_duplicates()
            df_na['type'] = 'intron'
            fdf = pd.concat([df, df_na]).fillna('-')
        fdf['coordinate'] = fdf.chr + ':' + fdf.coor.astype(str)
        fdf['exon_annotation'] = (fdf['strand'] + '|' + fdf['loc'] + '|' +
                                  fdf['gene'] + '|' + fdf['transcript'])
        fdf = fdf[['Fusion_id', 'coordinate', 'gene', 'type', 'exon_annotation']]
        fdf = fdf.sort_values(by='Fusion_id').reset_index(drop=True)
        return fdf

    start = time.process_time()
    fusion_bkpt_bedpe = utils.bed2bedpe(fusion_bkpt_df, 1, 0)
    bkpt5_df = fusion_bkpt_bedpe[['chr5', 'start5', 'end5', 'coor5', 'Fusion_id']]
    bkpt3_df = fusion_bkpt_bedpe[['chr3', 'start3', 'end3', 'coor3', 'Fusion_id']]

    exons_df = utils.read_df_or_file(config['exon_regions'], df=True)

    all_genes_df = utils.slop_bed(exons_df.copy(), 300)
    #all_genes_df.to_csv('exon_df.tsv', sep='\t', index=False)                          
    b5 = process_df(bkpt5_df, exons_df, all_genes_df)
    b3 = process_df(bkpt3_df, exons_df, all_genes_df)
    logger.info(f'exon annotation prep time: {round((time.process_time() - start)/60, 2)}mins')
    return b5, b3


class ProcessBreakpoints():
    """
    A class to annotate input fusion breakpoints

    Attributes
    ----------
    bkpt1, bkpt2 : str
        5' and 3' fusion breakpoint location eg:chr:765342
    exon_region_dict : dict
        Dict with exon start, end, body as keys and corresponding genomic positions
        as their values. Multiple transcripts of the same gene will be present.
    fusion_df : dataframe
       Df to store annotation of input fusion breakpoints

    Methods
    -------
    initate_processing(object):
        initate the annotation of fusion
    find_exon_loc(5' or 3' bkpt):
        find the part of exon the input bkpt belongs to and also the genomic region
    add2dict():
        format the data structure of the find_exon_loc()'s return'
    fusion_loc_annot_orient():
        fix bkpt orientation, and annotate based on exon location and bkpt genes
    fix_same_gene_fusion():
        find and annotate real and false same-gene fusions

    """

    def __init__(self, bkpt1, bkpt2, all_genes_df):
        """
        Parameters
        ----------
        bkpt1, bkpt2 : str
            5' and 3' fusion breakpoint location eg:chr:765342
        exon_region_dict : dict
            Dict with exon start, end, & body as keys and genomic pos as values
        """

        self.all_genes_df = all_genes_df
        self.b5_groupby = bkpt1.groupby('Fusion_id')
        self.b3_groupby = bkpt2.groupby('Fusion_id')
        self.COLS = ["Fusion_id", "5'-3'Gene_Partners", 'Fusion_Position',
                     'Fusion_Annotation', "5'_Exon_Annotation", "3'_Exon_Annotation",
                     "5'co-ordinate", "3'co-ordinate", 'Distance_between_breakpoints']

    def run_annotation(self, fusion_id):
        """
        Initates and manages the processing of fusion bkpts
        through other methods of the class.

        Return object
        -------------
        fusion_df : dataframe
            Stores the final annotations of input fusion in this object
        """
        start3 = time.process_time()
        self.bkpt_loc_dict = {'bkpt5': {}, 'bkpt3': {}}
        self.fusion_df = pd.DataFrame(columns=self.COLS)

        logger.debug(f'\n\n--------------Processing fusion_id: {fusion_id}---------------\n')
        logger.debug(f"5' breakpoint :\n")
        self.find_exon_loc(self.b5_groupby, fusion_id, direction='bkpt5')
        logger.debug(f"\n3' breakpoint :\n")
        self.find_exon_loc(self.b3_groupby, fusion_id, direction='bkpt3')

        #if len(self.bkpt_loc_dict) > 2:
            #logger.warning(f'Fusion is in 2+ exon locations{self.bkpt_loc_dict}')

        for count, fus_list in enumerate(self.fusion_loc_annot_orient()):
            logger.debug(f"\n\nfus_list:\n {fus_list}")
            fus_pos, bkpt1_list, bkpt2_list, fus_annot, gene_partners = fus_list

            # generates gene:transcript combinations
            fusion_df2 = pd.DataFrame(
                              list(itertools.product(bkpt1_list[0], bkpt2_list[0])),
                              columns=["5'_Exon_Annotation", "3'_Exon_Annotation"])
            # adding exon-position and coordinates to all trans-combinations
            fusion_df2 = fusion_df2.assign(Fusion_id=fusion_id,
                                           gene_partners=gene_partners,
                                           Fusion_Position=fus_pos,
                                           Fusion_Annotation=fus_annot,
                                           coordinate3=bkpt2_list[1],
                                           coordinate5=bkpt1_list[1],
                                           Distance_between_breakpoints=\
                                           self.get_genomic_distance(bkpt1_list[1],
                                                                     bkpt2_list[1]))
            # fusion in same genes needs secondary processing
            if fus_annot == 'same_gene':
                fusion_df2 = self.fix_same_gene_fusion(fusion_df2)

            # Renaming below column names by adding (')-prime here
            fusion_df2 = fusion_df2.rename(columns={
                                            'coordinate5': "5'co-ordinate",
                                            'coordinate3': "3'co-ordinate",
                                            'gene_partners': "5'-3'Gene_Partners"})

            self.fusion_df = pd.concat([self.fusion_df, fusion_df2], ignore_index=True)

            logger.debug(f'\nProcessed fusion_df (count:{count}): \n{self.fusion_df}\n')


    def find_exon_loc(self, bkpt_df, fusion_id, direction):
        """
        Fusion breakpoints(bkpt) occurs mostly in start/end/body of a gene's exon.
        This module determines the exon_loc and the genomic locations of input
        bkpt by searchching through each genomic position(vals) under each exon
        location (keys) of exon_region_dict. if not found in those three region,
        it assumes intronic region of the gene.

        Parameters
        ----------
        bkpt: str-list
            5' and 3' genomic location of fusion breakpoint  eg:['chr6', '765342']
        exon_region_dict : dict
            Dict with exon start, end, body as keys and corresponding genomic
            positions in bed format as dataframes
        out_df:
            output dataframe consists genomic positions, genes and transcript
            info in the .BED format
        bkpt_loc_dict: dictionary
            keys is the exon_location(start/stop/body/intron)
            value is a list with out_df and input_bkpt
        """
        fusion_id_df = get_group(bkpt_df, fusion_id)
        logger.debug(f'\nfusion_id: {fusion_id}, \nout_df: \n{fusion_id_df}\n\n')
        #if out_df.exon_annotation.notnull().any():
        try:
            for type_ in fusion_id_df.type.unique():
                logger.debug(f'Adding to bkpt_loc_dict - type_: {type_}\n')
                fusion_id_df_subset = fusion_id_df[fusion_id_df.type == type_]
                self.add_to_bkpt_loc_dict(''.join(type_),
                                          fusion_id_df_subset.exon_annotation, 
                                          fusion_id_df_subset.coordinate,
                                          direction)
        #else:
        except:
            # Not found in all 3 exon regions so adding empty series
            self.add_to_bkpt_loc_dict('intron', pd.Series('-|NA|-|-'),
                                      fusion_id_df.coordinate, direction)

    def add_to_bkpt_loc_dict(self, loc, loc_df, bkpt, direction):
        bkpt = ''.join(bkpt.unique())

        if loc in self.bkpt_loc_dict[direction]:
            self.bkpt_loc_dict[direction][loc + '2'] = [loc_df, bkpt]
        else:
            self.bkpt_loc_dict[direction][loc] = [loc_df, bkpt]

        logger.debug(f'\nUpdated bkpt_loc_dict: \n{self.bkpt_loc_dict}\n')


        """
        bkpt_loc_dict {'end':
                                [
                                 0       +|End_E6|EML4|NM_001145076
                                 1          +|End_E7|EML4|NM_019063
                                 2    +|End_E8|EML4|ENST00000401738
                                 Name: exon_annotation, dtype: object, 'chr2:42508112'
                                 ],
                       'body':
                                [
                                 0    -|Body_E20|ALK|NM_004304
                                 Name: exon_annotation, dtype: object, 'chr2:29446376'
                                 ]
                    }

        """

    @staticmethod
    def get_genomic_distance(bkpt1_coord, bkpt2_coord):
        distance = 'NA'
        bkpt1_list = bkpt1_coord.split(':')
        bkpt2_list = bkpt2_coord.split(':')
        # if both bkpt in same chromosome, calculate distance b/w them
        if bkpt1_list[0] == bkpt2_list[0]:
            distance = int(max(int(bkpt1_list[1]), int(bkpt2_list[1])) -
                           min(int(bkpt1_list[1]), int(bkpt2_list[1])))
        return distance

    def fusion_loc_annot_orient(self):
        """
        Predicts 'fus_annotation' column based on the exon_locations of both bkpts

        The breakpoint reported by fusion callers may not always be in correct
        5'-3' orientations. Depending on the exon_locs predicted from prev step,
        some biologically applicable annotation columns like "fusion-location" and
        "fusion-annotation" column values and correct orientation are determined and
        applied while passing values to check_same_gene(). The check_same_gene()
        checks if fusion is happening in the same gene, if so, the fusion-annotation
        is modified to 'same-gene' as they need futher processing in next the step.
        It also adds "5'-3' gene partners" to the output.

        Parameters & functions
        ----------------------
        bkpt: str-list
            5' and 3' genomic location of fusion breakpoint  eg:['chr6', '765342']
        fus_annot : str
            predicted fusion-annot term based on exon locations
        exon_annots:
            exon_level annotations like genomic positions, strand, genes and
            transcript identified by the prev method find_exon_loc()
        bkpt_loc_dict/bdict: dictionary
            keys is the exon_location(start/stop/body/intron)
            value is a list with out_df(from prev step) and bkpt-coordinates
        check_same_gene(str, list, list):
            check if fusion are in the same chromosome and gene
        """

        def check_same_gene(fus_pos, fus_annot, bkpt1_list, bkpt2_list):
            """
            Generator function that checks if the input fusion is
            in the same gene. If there are two genes in the same
            genomic coordinate of the fusion breakpoint, this generator,
            yeilds each one of them
            """

            fus_annot_list = []
            exon5_annot = bkpt1_list[0].str.split('|').str
            exon3_annot = bkpt2_list[0].str.split('|').str
            # a copy to preserve overwriting
            bkpt1_series = bkpt1_list[0]
            bkpt2_series = bkpt2_list[0]
            # looping to add multiple genes at same genomic co-ordinates"
            for gene5 in exon5_annot[2].unique():
                for gene3 in exon3_annot[2].unique():

                    bkpt1_list[0] = self.get_gene_annot(bkpt1_series, gene5)
                    bkpt2_list[0] = self.get_gene_annot(bkpt2_series, gene3)
                    gene_partners = gene5 + '-' + gene3
                    #gene_partners =  f'<a href="{config["genecards_url"]}{gene5}">{gene5}</a>-<a href="{config["genecards_url"]}{gene3}">{gene3}</a>'

                    # same chr and same gene
                    if gene5 == gene3 and gene5 == gene3:
                        fus_annot = 'same_gene'

                    logger.debug(f"\n\ncheck_same_gene: \n\n fus_pos:\n\n'{fus_pos}'\n, bkpt1_list:\n{bkpt1_list}\n, bkpt2_list:\n\n{bkpt2_list}\n, fus_annot:\n\n{fus_annot}\n, gene_partners:\n\n{gene_partners}\n")
                    fus_annot_list.append([fus_pos, bkpt1_list, bkpt2_list,
                                            fus_annot, gene_partners])
            
            return fus_annot_list

        def determine_fus_annot_and_pos(bdict):

                # bkpt at exon start and end is always at exon end2start orientation
                if all(loc in bdict for loc in ['end', 'start']):
                    fus_annot = 'Fusion-Candidate'
                    return check_same_gene('Exon-Exon_boundary', fus_annot, bdict['end'], bdict['start'])
                # bkpt at exon body and end is always at exon end2body orientation
                elif all(loc in bdict for loc in ['end', 'body']):
                    fus_annot = 'Fusion-Candidate'
                    return check_same_gene('Exon_boundary-Exon_Body', fus_annot, bdict['end'], bdict['body'])
                # bkpt at exon body and start is always body2start orientation
                elif all(loc in bdict for loc in ['body', 'start']):
                    fus_annot = 'Fusion-Candidate'
                    return check_same_gene('Exon_Body-Exon_boundary', fus_annot, bdict['body'], bdict['start'])

                elif all(loc in bdict for loc in ['start', 'start2']):
                    fus_annot = 'Fused_Start-to-Start'
                    return check_same_gene('Exon-Exon_boundary', fus_annot, bdict['start'], bdict['start2'])

                elif all(loc in bdict for loc in ['end', 'end2']):
                    fus_annot = 'Fused_End-to-End'
                    return check_same_gene('Exon-Exon_boundary', fus_annot, bdict['end'], bdict['end2'])

                elif all(loc in bdict for loc in ['body', 'body2']):
                    fus_annot = 'Fused_Body-to-Body'
                    return check_same_gene('Exon_Body-Exon_Body', fus_annot, bdict['body'], bdict['body2'])
                # bkpt at intron and exon-start is always body2start orientation
                elif all(loc in bdict for loc in ['intron', 'start']):
                    fus_annot = 'Potential_Fusion-Candidate'
                    return check_same_gene('Intron-Exon_boundary', fus_annot, bdict['intron'], bdict['start'])
                # bkpt at intron and exon-end is always end2intron orientation
                elif all(loc in bdict for loc in ['intron', 'end']):
                    fus_annot = 'Potential_Fusion-Candidate'
                    return check_same_gene('Exon_boundary-Intron', fus_annot, bdict['end'], bdict['intron'])
                # bkpt at intron and exon-body is always intron2end orientation
                elif all(loc in bdict for loc in ['intron', 'body']):
                    fus_annot = 'Potential_Fusion-Candidate'
                    return check_same_gene('Intron-Exon_Body', fus_annot, bdict['intron'], bdict['body'])   
                # both bkpt annotated as intron not found found in exons
                elif all(loc in bdict for loc in ['intron', 'intron2']):
                    fus_annot = 'Fusion-Candidate'
                    return check_same_gene('within_Gene_but_not_Exon_boundary', fus_annot, bdict['intron'], bdict['intron2'])
                else:
                    logger.warning(f'trap- fusion doesn\'t fall into any exon region -{bdict}')

        for bdict5_loc in self.bkpt_loc_dict['bkpt5']:
            for bdict3_loc in self.bkpt_loc_dict['bkpt3']:

                bdict = {bdict5_loc: self.bkpt_loc_dict["bkpt5"][bdict5_loc]}
                threep_data = self.bkpt_loc_dict["bkpt3"][bdict3_loc]

                # handling same exon locations for both bkpts
                bdict3_loc = bdict3_loc + '2' if bdict3_loc in bdict else bdict3_loc
                bdict[bdict3_loc] = threep_data
                yield from determine_fus_annot_and_pos(bdict)



    @staticmethod
    def get_gene_annot(bkpt_series, gene):
        ''' Splitting bkpt series on gene '''
        return bkpt_series[bkpt_series.str.contains('|' + gene + '|', regex=False)]

    @staticmethod
    def fix_same_gene_fusion(fusion_df):
        """
        Further process fusion in same-gene to determine true fusion event.

        Fusion reported as 'same-gene' from prev step is not always true fusion.
        This method detrmine true fusion by checking the gene-transcripts and
        exon numbers and annotates them accordingly.

        Parameters & fuctions
        ---------------------
        fusion_df: Dataframe
            contains annotations from prev steps for the fusion bkpts.
        fusion_conditions(Series, Series) :
            applies multiple condition on each line from input fusion_df to
            determine and annotate same-gene fusion accordingly
        exon_annot:
            exon_level annotation column from fusion_df
        """
        exon5_annot = fusion_df["5'_Exon_Annotation"].str.split('|')
        exon3_annot = fusion_df["3'_Exon_Annotation"].str.split('|')

        def fusion_conditions(exon5_annot, exon3_annot):
            # check if fusions are from same transcripts and neighbor exons
            same_trans = exon5_annot[3] == exon3_annot[3]
            exon_calc = abs(int(exon3_annot[1].split('E')[-1]) -
                         int(exon5_annot[1].split('E')[-1])) 
            #print(f"\n fusion_conditions: {exon3_annot[1]}, {exon5_annot[1]}, {exon_calc} ")

            if same_trans:

                if exon_calc == 0:
                    # same gene same exon
                    return 'Within_Gene_same_Exons'
                elif exon_calc == 1:
                    # if immidiate neighbor exon
                    return 'Within_Gene_Neighboring_Exons'
                elif exon_calc > 1:                     
                    # same gene any exon but same/neighbor exon
                    return 'Same_Gene_Fusion'
            
            else:
                return 'remove'
        try:
            func = np.vectorize(fusion_conditions)
            fusion_df["Fusion_Annotation"] = func(exon5_annot, exon3_annot)
        except ValueError:
            print("****Error in fusion_conditions***")
            fusion_df["Fusion_Annotation"] = 'remove'

        return fusion_df[fusion_df.Fusion_Annotation != 'remove']


def get_igv_link(row):
    ''' Adding IGV URL for each breakpoint '''

    def get_val(annot, coordinate):
        if annot == 'NA' or annot == 'nan':
            return coordinate
        annot_split = annot.split('|')
        final_annot = annot_split[2] + ':' + annot_split[3] + \
                       ':' + annot_split[1].split('_')[1]
        return str(final_annot)

    five_val = get_val(row["5'_Exon_Annotation"], row["5'co-ordinate"])
    three_val = get_val(row["3'_Exon_Annotation"], row["3'co-ordinate"])
    igv_id = five_val + "-" + three_val
    igv_link = config['igv_link'].replace("IGV_ID", igv_id.strip())
    return igv_link


if __name__ == '__main__':
    main()
