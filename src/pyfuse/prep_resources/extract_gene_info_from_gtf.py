#!/usr/bin/env python3


import os
import sys
import hashlib
import argparse
import logging
import validators
import numpy as np
import pandas as pd
import urllib.request
import re
from tqdm import tqdm
from pathlib import Path

from urllib.parse import urlparse, urljoin
from pyfuse.utils.common_utils import utils
tqdm.pandas()

config = utils.get_config_vars()

resource_files_dict = {
    "exon_regions": "genes_by_location.tsv.gz",
    "all_gene_exons": "all_genes_exons_transcripts.bed.gz",
    "start_codon_file": "start_codons.txt.gz",
}

DEFAULT_GTF_URL_GRCH37 = config['default_refseq_gtf_url_grch37']
DEFAULT_ANNOT_SUMMARY_GRCH37 = config['default_refseq_assembly_summary_url_grch37']
DEFAULT_GTF_URL_GRCH38 = config['default_refseq_gtf_url_grch38']
DEFAULT_ANNOT_SUMMARY_GRCH38 = config['default_refseq_assembly_summary_url_grch38']
DEFAULT_MANE_URL = config['default_mane_url']

logger = logging.getLogger(__name__)


def display_genome(genome: str) -> str:
    mapping = {
        'grch37': 'GRCh37',
        'grch38': 'GRCh38',
    }
    return mapping.get(genome, genome)


class RequiredHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Append a required marker to required optional arguments in help output."""

    def _get_help_string(self, action):
        help_text = super()._get_help_string(action)
        if action.required and action.option_strings and '[required]' not in help_text.lower():
            return f"[required] {help_text}"
        return help_text


def add_arguments(parser):
    parser.add_argument('--genome', dest='genome', required=True,
                        choices=['grch37', 'grch38'],
                        help='reference genome build used by input annotations')
    parser.add_argument('-s', dest='source', required=False, help='source of the gene annotation file',
                        choices=['refseq-gtf', 'ucsc-gtf'],
                        default='refseq-gtf')
    parser.add_argument('-g', dest='gtf', required=False,
                        help='reference gtf file (auto-switches to GRCh38 default if --genome grch38 and -g not provided)',
                        default=None)
    parser.add_argument('-a', dest='annot_summary', required=False,
                        help='assembly summary report for refseq (auto-switches to GRCh38 default if --genome grch38 and -a not provided)',
                        default=None)
    parser.add_argument('-o', dest='out_path', help='output path', required=True, type=(os.path.abspath))
    parser.add_argument('-n', dest='out_name', help='output folder name', required=False, type=(os.path.abspath))
    parser.add_argument('-t', dest='roi', type=(os.path.abspath), required=False, help='target roi bed to filter annotations')
    parser.add_argument('-d', dest='annot_delim', default=' ', help='delimiter in annotation column of gtf(last col)')
    parser.add_argument('-l', dest='long_trans', action='store_true', help='flag to only choose longest transcript per gene')
    parser.add_argument('-c', dest='req_annot_cols', help='required annotation column of gtf(last col)',
                        default=['transcript_id', 'gene_id', 'gene', 'exon_number', 'exonCount'])
    parser.add_argument('-m', '--mane-file', '--mane-url', dest='mane_file', required=False,
                        default=DEFAULT_MANE_URL,
                        help='MANE summary local file path or URL; URL can be a MANE directory or direct summary file')
    return parser


def warn_ucsc_source_usage() -> None:
    logger.warning(
        "You selected source='ucsc-gtf'. "
        "UCSC permits broad reuse for most core data, but some tracks/datasets have "
        "separate restrictions from upstream providers. "
        "Before redistribution, verify dataset-specific terms and avoid restricted tracks."
    )


def parse_args():
    parser = argparse.ArgumentParser(description='Script 2', formatter_class=RequiredHelpFormatter)
    add_arguments(parser)
    return parser.parse_args()

def mkdir(path):
    os.makedirs(path, exist_ok=True)


def intersect_gtf_to_roi(gtf_df: pd.DataFrame, roi_bed_path: str) -> pd.DataFrame:
    """Keep annotation rows that overlap any interval in a ROI BED file."""
    roi_df = pd.read_csv(
        roi_bed_path,
        sep='\t',
        comment='#',
        header=None,
        usecols=[0, 1, 2],
        names=['chrom', 'start', 'end'],
        dtype={0: str},
    )

    if roi_df.empty:
        logger.warning("ROI BED is empty; returning unfiltered annotations")
        return gtf_df

    roi_df['start'] = pd.to_numeric(roi_df['start'], errors='coerce')
    roi_df['end'] = pd.to_numeric(roi_df['end'], errors='coerce')
    roi_df.dropna(subset=['start', 'end'], inplace=True)
    roi_df['start'] = roi_df['start'].astype(int)
    roi_df['end'] = roi_df['end'].astype(int)

    a = gtf_df.copy()
    a['start'] = pd.to_numeric(a['start'], errors='coerce').astype('Int64')
    a['end'] = pd.to_numeric(a['end'], errors='coerce').astype('Int64')
    a.dropna(subset=['start', 'end'], inplace=True)
    a['start'] = a['start'].astype(int)
    a['end'] = a['end'].astype(int)

    if a.empty:
        return a

    a_cols = list(a.columns)
    a_bed = utils.df2pbt(a)
    roi_bed = utils.df2pbt(roi_df[['chrom', 'start', 'end']])
    overlapped = a_bed.intersect(roi_bed, u=True).to_dataframe(
        header=None,
        disable_auto_names=True,
    )

    if overlapped.empty:
        logger.warning("No annotations overlapped the provided ROI BED; returning empty table")
        return a.iloc[0:0].copy()

    overlapped.columns = a_cols
    return overlapped


def normalize_tx_id(tx):
    tx = str(tx).strip()
    if not tx or tx.lower() == 'nan':
        return ''
    return tx.split('.')[0]


def resolve_local_or_url(path_or_url, file_path, label):
    if not path_or_url:
        return None
    if os.path.exists(path_or_url):
        logger.info("Using local %s: %s", label, path_or_url)
        return path_or_url
    if validators.url(path_or_url):
        local_file = os.path.join(file_path, os.path.basename(urlparse(path_or_url).path))
        if not local_file or local_file.endswith('/'):
            sys.exit(f"\n[ERROR] {label} URL should resolve to a file, got directory URL: {path_or_url}")
        if not os.path.exists(local_file):
            logger.info("Downloading %s from %s", label, path_or_url)
            urllib.request.urlretrieve(path_or_url, local_file)
        else:
            logger.info("Reusing cached %s: %s", label, local_file)
        return local_file
    sys.exit(f"\n[ERROR] {label} input {path_or_url} is neither a valid file nor a valid url")


def parse_version(version_txt):
    nums = re.findall(r'\d+', version_txt)
    return tuple(int(x) for x in nums)


def resolve_mane_summary_url(mane_url):
    if mane_url.endswith('.summary.txt.gz'):
        return mane_url

    listing_url = mane_url if mane_url.endswith('/') else mane_url + '/'
    try:
        with urllib.request.urlopen(listing_url) as resp:
            html = resp.read().decode('utf-8', errors='ignore')
    except Exception as e:
        sys.exit(f"\nERROR: Unable to read MANE directory URL {listing_url}: {e}")

    matches = re.findall(r'(MANE\.GRCh38\.v([0-9\.]+)\.summary\.txt\.gz)', html)
    if not matches:
        sys.exit(f"\nERROR: No MANE GRCh38 summary file found at {listing_url}")

    best_name, _ = max(matches, key=lambda x: parse_version(x[1]))
    return urljoin(listing_url, best_name)


def confirm_user_choice(prompt_message):
    while True:
        try:
            ans = input(prompt_message).strip().lower()
        except EOFError:
            return False
        if ans in {'yes', 'y'}:
            return True
        if ans in {'no', 'n'}:
            return False
        print("Please answer yes or no.")


def get_default_refseq_urls(genome):
    if genome == 'grch38':
        return DEFAULT_GTF_URL_GRCH38, DEFAULT_ANNOT_SUMMARY_GRCH38
    return DEFAULT_GTF_URL_GRCH37, DEFAULT_ANNOT_SUMMARY_GRCH37


def load_mane_transcript_map(mane_file):
    mane_df = pd.read_csv(mane_file, sep='\t', dtype=str, compression='infer')

    if 'MANE_status' not in mane_df.columns:
        sys.exit("\nERROR: MANE summary file missing required column 'MANE_status'")

    gene_candidates = ['symbol', 'HGNC_symbol', 'Gene', 'gene', 'Gene name']
    gene_col = next((c for c in gene_candidates if c in mane_df.columns), None)
    name_candidates = ['name', 'Name', 'description', 'Description']
    name_col = next((c for c in name_candidates if c in mane_df.columns), None)

    tx_candidates = ['RefSeq_nuc', 'RefSeq_nuc_acc']
    tx_col = next((c for c in tx_candidates if c in mane_df.columns), None)
    if not tx_col:
        sys.exit(
            "\nERROR: MANE summary file missing RefSeq transcript column. Expected one of: "
            "RefSeq_nuc, RefSeq_nuc_acc"
        )

    classes = sorted([x for x in mane_df['MANE_status'].dropna().astype(str).unique().tolist() if x.strip()])
    if not classes:
        sys.exit("\nERROR: No MANE classes found in MANE_status column")

    subset_cols = ['MANE_status', tx_col] + ([gene_col] if gene_col else []) + ([name_col] if name_col else [])
    sub = mane_df[subset_cols].copy()
    sub = sub.rename(columns={tx_col: 'transcript_id'})
    sub['transcript_id'] = sub['transcript_id'].astype(str).map(normalize_tx_id)
    sub = sub[(sub['transcript_id'] != '') & (sub['MANE_status'].notna())]
    sub['gene'] = sub[gene_col].astype(str) if gene_col else '.'
    sub['gene'] = sub['gene'].replace({'nan': '.'})
    sub['gene_name'] = sub[name_col].astype(str) if name_col else '.'
    sub['gene_name'] = sub['gene_name'].replace({'nan': '.'})

    mane_long = sub[['gene', 'transcript_id', 'gene_name', 'MANE_status']].drop_duplicates()
    tx_info = (mane_long[['transcript_id', 'gene', 'gene_name']]
               .drop_duplicates()
               .sort_values(by=['gene', 'transcript_id', 'gene_name'])
               .groupby('transcript_id', as_index=False)
               .first())

    matrix_df = (mane_long.assign(value=1)
                 .pivot_table(index=['transcript_id'],
                              columns='MANE_status',
                              values='value',
                              aggfunc='max',
                              fill_value=0)
                 .reset_index())

    matrix_df = tx_info.merge(matrix_df, on='transcript_id', how='left')

    for class_name in classes:
        if class_name not in matrix_df.columns:
            matrix_df[class_name] = 0

    matrix_df = matrix_df[['gene', 'gene_name', 'transcript_id'] + classes]
    matrix_df[classes] = matrix_df[classes].astype(int)
    matrix_df = matrix_df.sort_values(by=['gene', 'transcript_id']).reset_index(drop=True)

    matrix_df = matrix_df.rename(columns={
        'gene': 'gene_symbol',
        'gene_name': 'gene_function',
        'transcript_id': 'transcript'
    })

    metadata = {
        'gene_column': gene_col if gene_col else 'NA',
        'gene_name_column': name_col if name_col else 'NA',
        'transcript_column': tx_col,
        'mane_classes': ','.join(classes),
        'row_count_input': str(len(mane_df)),
        'row_count_matrix': str(len(matrix_df))
    }
    return matrix_df, metadata


def main(args):
    if not logging.getLogger().handlers:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s [%(module)s]: %(message)s")

    genome_label = display_genome(args.genome)

    if args.source == 'ucsc-gtf':
        warn_ucsc_source_usage()

    default_gtf_url, default_annot_summary_url = get_default_refseq_urls(args.genome)
    gtf_input = args.gtf if args.gtf else default_gtf_url
    annot_summary_input = args.annot_summary if args.annot_summary else default_annot_summary_url

    logger.info("Preparing PyFuse resources for genome '%s' using source '%s'", genome_label, args.source)

    if args.gtf:
        logger.info("Using user-provided GTF source: %s", args.gtf)
    else:
        logger.info("Using default GTF source from settings: %s", default_gtf_url)

    if args.annot_summary:
        logger.info("Using user-provided annotation summary source: %s", args.annot_summary)
    else:
        logger.info("Using default annotation summary source from settings: %s", default_annot_summary_url)

    parent_out_path = os.path.join(args.out_path, 'custom_resource_files')
    mkdir(parent_out_path)


    if not args.out_name:
        out_path = os.path.join(parent_out_path, os.path.basename(gtf_input).split('.')[0])
    else:
        out_path = os.path.join(parent_out_path, args.out_name)
    
    mkdir(out_path)

    required_chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
                       'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                       'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                       'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']
    file_path = os.path.join(out_path, 'downloaded_gtf_files')
    mkdir(file_path)

    if os.path.exists(gtf_input):
        logger.info("Using local GTF file: %s", gtf_input)
        gtf_file = gtf_input
        annot_summary_file = annot_summary_input
        if annot_summary_input:
            if os.path.exists(annot_summary_input):
                logger.info("Using local annotation summary file: %s", annot_summary_input)

    elif validators.url(gtf_input):
        gtf_name = os.path.basename(urlparse(gtf_input).path)
        summary_name = os.path.basename(urlparse(annot_summary_input).path)
        gtf_file = os.path.join(file_path, gtf_name)
        annot_summary_file = os.path.join(file_path, summary_name)

        if os.path.exists(gtf_file) and os.path.exists(annot_summary_file):
            logger.info("Reusing cached downloads for GTF and annotation summary from %s", file_path)
        else:
            logger.info("Downloading reference inputs into %s", file_path)
            urllib.request.urlretrieve(gtf_input, gtf_file)
            urllib.request.urlretrieve(annot_summary_input, annot_summary_file)
    else:
        sys.exit(
            "\n[ERROR] GTF/annotation-summary inputs are invalid. "
            f"Received: gtf={gtf_input}, annot_summary={annot_summary_input}. "
            "Provide valid local file paths or valid URLs."
        )

    #output_file = os.path.join(out_path, os.path.basename(gtf_file).split('.')[0] + '_formatted.tsv')


    # Reading  gtf
    if args.source == 'refseq-gtf':
        if not annot_summary_input:
            sys.exit(
                "\n[ERROR] RefSeq GTF parsing requires an annotation summary report (-a/--annot-summary). "
                "It is usually available in the same source directory as the GTF."
            )
        annot_summary_file = resolve_local_or_url(annot_summary_file, file_path, label='annotation summary')

    logger.info("Reading GTF annotations")
    raw_refseq_df = pd.read_csv(gtf_file, sep="\t", comment='#',
                                usecols=[0, 2, 3, 4, 6, 8], header=None)
    raw_refseq_df = raw_refseq_df[raw_refseq_df.iloc[:, 1].isin(['exon', 'transcript', 'start_codon'])]
    raw_refseq_df.reset_index(drop=True, inplace=True)
    logger.info("Retained %s rows after selecting exons, transcripts, and start codons", len(raw_refseq_df))

    # Process each annotation and extract the required columns and merge to main df
    def process_column(df, column):
        df['key'] = df[column].str.strip().str.split(args.annot_delim).str[0]
        df['val'] = df[column].str.strip().str.split(args.annot_delim).str[-1].str.replace('"', "")
        df = df.pivot(columns='key', values='val')
        #print(df)
        return df

    annots_df = pd.DataFrame()
    annots_cols = raw_refseq_df[8].str.split(";", expand=True)
    raw_refseq_df.drop(raw_refseq_df.columns[-1], axis=1, inplace=True)
    req_cols = [ "chrom", "region", "start", "end", "strand", 'ID',
                 "gene_id", "gene", "transcript", "exonCount", 'exon_number']
    #chrom_list = ["chr" + str(i) for i in list(range(1, 23)) + ["X", "Y"]]
    for col in tqdm(annots_cols.columns):
        col_result = process_column(annots_cols, col)
        for col in col_result.columns:
            #print('\neachcol of colresult', col)
            if col in args.req_annot_cols:
                if col not in annots_df:
                    annots_df[col] = col_result[col]
                else:# Merge on the index to ensure that rows are aligned correctly
                    col_result.dropna(subset=[col], inplace=True)
                    annots_df[col] = annots_df[col].combine_first(col_result[col])
                #print('\n\nIN:  ',col, '\n', col_result)
    annots_df.rename(inplace=True, columns={'gene_id': 'gene',
                                            'exonStarts': 'start',
                                            'exonEnds': 'end',
                                            'transcript_id': 'transcript',
                                            'exonCount': 'exon_number'})
    annots_df = drop_dup_col(annots_df)
    # print(annots_df, raw_refseq_df)


    raw_refseq_df = pd.concat([raw_refseq_df, annots_df], axis=1)
    req_cols = req_cols[:5] + [i for i in req_cols[5:] if i in raw_refseq_df]
    raw_refseq_df.columns = req_cols

    # find transcript rows, extract and cmove their coordinates as separate columns
    # only include the long transcript per gene based on user input in -l
    def reasign_transcripts(gene_df):
        transdf = gene_df.loc[gene_df['region'] == 'transcript'] \
                             [['transcript_id', 'start', 'end']]
        gene_df.drop(transdf.index, inplace=True)
        transdf["trans_len"] = transdf["end"] - transdf["start"]
        transdf.reset_index(inplace=True, drop=True)
        if len(transdf) > 1:
            transdf = transdf.iloc[transdf["trans_len"].idxmax()]
            gene_df = gene_df.loc[gene_df["transcript_id"] == transdf["transcript_id"]]
        gene_df[['txStart', 'txEnd']] = int(transdf["start"]), int(transdf["end"])
        gene_df = gene_df[['chrom', 'start',  'end', 'strand', 'exonCount',
                           'gene_id', 'transcript_id', 'txStart','txEnd']]
        return gene_df

    logger.info("Preparing transcript-level annotations")
    if args.long_trans:
        final_df = raw_refseq_df.groupby(["gene_id"], group_keys=False).progress_apply(lambda x: reasign_transcripts(x))
    else:
        final_df = raw_refseq_df

    # Convert refseq accession id to ucsc style chrom number
    if args.source == 'refseq-gtf':
        logger.info("Converting RefSeq accession identifiers to UCSC chromosome names")
        with open(annot_summary_file, 'r') as fh:
            for header_index, l in enumerate(fh):
                if 'UCSC-style-name' in l:
                    break

        cols = final_df.columns
        assembly_df = pd.read_csv(annot_summary_file, sep="\t", header=header_index,
                                  usecols=['RefSeq-Accn', 'UCSC-style-name'])
        final_df = final_df.merge(assembly_df, left_on='chrom', right_on='RefSeq-Accn', how='left') \
                           .drop(['chrom','RefSeq-Accn'], axis=1)\
                           .rename(columns={'UCSC-style-name': 'chrom'})
        final_df.dropna(subset=['chrom'], inplace=True)
        final_df = final_df[final_df['chrom'].isin(required_chroms)]
        final_df.reset_index(inplace=True)
        #print(final_df)
        final_df = final_df[cols]

        logger.info("Retained %s rows after UCSC chromosome filtering", len(final_df))

    if args.roi:
        logger.info("Filtering annotations to the provided ROI BED regions")
        final_df = intersect_gtf_to_roi(final_df, args.roi)

    # remove unassigned transcripts
    final_df = final_df[~final_df['transcript'].str.contains('unassigned', na=False)]
    #final_df.to_csv(output_file, sep='\t', index=False)
    #print(final_df)
    #print('\nOUTPUT: ', output_file)


    ####################################################
    #   all_genes_exons_transcripts.bed
    ####################################################
    logger.info("Writing exon-level annotation resources")
    # Columns to concatenate
    all_genes_cols = ['chrom', 'start', 'end', 'gene', 'annotation']
    cols_to_concat = ['strand', 'exon_number', 'gene', 'transcript']

    def concat_with_pipe(row):
        return '|'.join(str(x) for x in row)
    # Concatenate and create a new column
    all_genes_df = final_df.loc[final_df["region"] == "exon"].copy()
    all_genes_df['exon_number'] = 'E' +  all_genes_df['exon_number'].astype(str)

    annotations = []
    for row in tqdm(all_genes_df[cols].itertuples(index=False, name=None),
                    total=len(all_genes_df),
                    desc="Building annotation"):
        annotations.append("|".join("" if x is None else str(x) for x in row))

    all_genes_df["annotation"] = annotations


    all_genes_df = all_genes_df[all_genes_cols]
    all_genes_df.dropna(inplace=True)
    all_genes_df.reset_index(drop=True, inplace=True)

    #resource_files_dict['all_gene_exons'] = resource_files_dict['all_gene_exons'] + '.zip'
    write_table(all_genes_df, out_dir=out_path, filename=resource_files_dict['all_gene_exons'],
                sep="\t", index=False)



    ####################################################
    #  Code for start codon file:
    ####################################################
    
    def find_start_codon(gene_df):
        # Initialize the new column with zeros
        gene_df.insert(loc=1, column='start_codon', value='-')
        if gene_df['region'].str.contains('start_codon').any():
            # get start codon start pos by strand for first cds found fo rthe transcript
            st_cdn_row = gene_df.loc[gene_df['region'] == 'start_codon'].iloc[0][['strand', 'start', 'end']]
            start_codon_start = st_cdn_row['start'] if st_cdn_row['strand'] == '+' else st_cdn_row['end']
            # Identify the index of the exon row where the CDS start fits within the range
            start_codon_exon_idx = gene_df.index[(gene_df['region'] == 'exon')
                                                & (gene_df['start'] <= start_codon_start)\
                                                & (gene_df['end'] >= start_codon_start)].tolist()[0]
            # If a hit was found, set the value from that row onwards to 1

            if start_codon_exon_idx is not None:
                gene_df.loc[(gene_df.index >= start_codon_exon_idx) & (gene_df['region'] == 'exon'), 'start_codon'] = str(start_codon_start)

        gene_df = gene_df[gene_df['region'] == 'exon'].drop(columns=['region'])
        exon_count_col = gene_df.pop('exon_number')
        gene_df.insert(loc=4, column='length', value=(abs(gene_df['start'] - gene_df['end'])+1))
        gene_df.insert(6, column='exon_number', value='E' + exon_count_col.astype(str))
        return gene_df

    logger.info("Determining start codon positions")
    # = final_df.groupby(['gene', 'transcript'], group_keys=False).apply(lambda x: find_start_codon(x))
    #start_codon_df.reset_index(drop=True, inplace=True)
    
    g = final_df.groupby(['gene', 'transcript'], group_keys=False, sort=False)

    start_codon_df = pd.concat(
        [find_start_codon(x) for _, x in tqdm(g, total=g.ngroups, desc="Start codon")],
        ignore_index=True
    )

    write_table(start_codon_df, out_dir=out_path, filename=resource_files_dict['start_codon_file'],
                sep="\t", index=False)
 

    ####################################################
    #   genes_by_location.tsv
    ####################################################
    logger.info("Writing exon-location resource")
    genes_by_loc_padding = 15

    def split_by_loc(row):
        row = row.copy()
        start = row['start']
        end = row['end']
        strand = row['strand']
        row.drop(['start', 'end'], inplace=True)
        pad = 30 
        # Calculate the original region length
        exon_len = end - start + 1
        start_type, end_type = ('End', 'Start') if strand == '-' else ('Start', 'End')
  
        if exon_len < 32:
            # if exon len is below 32 return exon as body cuz it doesn't make snse to split any further
            new_df = pd.DataFrame([{'start': start, 'end': end, 'type': 'body'}])
        
        else:
            # Adjust padding if the original region length is less than 60
            if 32 <= exon_len <= 65 :
                pad = 15
            # Create three new rows based on the described logic
            row_start = {'start': start, 'end': start + (pad-1), 'type': start_type}
            row_body = {'start': start + pad, 'end': end - pad, 'type': 'Body'}
            row_end = {'start': end - (pad-1), 'end': end, 'type': end_type}\
            

            # Create a DataFrame from the new rows
            new_df = pd.DataFrame([row_start, row_body, row_end])

        row_reps = pd.DataFrame(np.tile(row.values, len(new_df)).reshape(len(new_df), -1),
                                    columns=row.index)
        # Return a DataFrame
        return pd.concat([new_df, row_reps], axis=1)


    # Apply the function to split the rows and reconstruct the DataFrame
    genes_by_loc_df = final_df.loc[final_df['region'] == 'exon'].copy()
    genes_by_loc_df['start'] = genes_by_loc_df['start'] - genes_by_loc_padding
    genes_by_loc_df['end'] = genes_by_loc_df['end'] + genes_by_loc_padding

    result_list = [
        split_by_loc(row)
        for _, row in tqdm(genes_by_loc_df.iterrows(), total=len(genes_by_loc_df), desc="gene_by_exon_location")
    ]
    genes_by_loc_df = pd.concat(result_list, ignore_index=True)

    # adding prefix 'E' to exon number
    genes_by_loc_df['loc'] = genes_by_loc_df['type'].astype(str) + '_E' + genes_by_loc_df['exon_number'].astype(str)
    genes_by_loc_df.drop(columns=['exon_number'], inplace=True)
    genes_by_loc_df['type'] = genes_by_loc_df['type'].astype(str).str.lower()
    genes_by_loc_cols = ['chrom', 'start', 'end', 'strand', 'loc',
                         'gene', 'transcript', 'type']
    genes_by_loc_df = genes_by_loc_df[genes_by_loc_cols]
    genes_by_loc_df.reset_index(drop=True, inplace=True)
    
    write_table(genes_by_loc_df, out_dir=out_path, filename=resource_files_dict['exon_regions'],
                sep="\t", index=False)

    mane_enabled = False
    mane_input = args.mane_file or DEFAULT_MANE_URL
    mane_summary_source, mane_file = None, None
    mane_meta = {}

    if mane_input == DEFAULT_MANE_URL:
        if args.genome == 'grch38':
            mane_enabled = True
        else:
            logger.info("Skipping MANE resource generation for %s", genome_label)
    else:
        if args.genome == 'grch37':
            warn_msg = (
                "\nWARNING: MANE is primarily intended for GRCh38. "
                "You are running grch37 with a user-provided MANE input. "
                "If this is a back-mapped grch37 resource, annotations may be incomplete or incorrect. "
                "Proceed? [yes/no]: "
            )
            if validators.url(mane_input) and re.search(r'grch37|hg19', mane_input, flags=re.IGNORECASE):
                warn_msg = (
                    "\nWARNING: You provided a GRCh37/HG19 MANE URL. Official MANE is GRCh38-centered; "
                    "back-mapped GRCh37 resources may have missing or incorrect annotations. "
                    "Are you sure this is the intended back-mapped version and do you want to proceed? [yes/no]: "
                )
            if not confirm_user_choice(warn_msg):
                logger.info("Skipping MANE resource generation based on user response")
                mane_enabled = False
            else:
                mane_enabled = True
        else:
            mane_enabled = True

    if mane_enabled:
        logger.info("Writing MANE annotation resource")
        if os.path.exists(mane_input):
            mane_summary_source = mane_input
            logger.info("Using local MANE summary file: %s", mane_summary_source)
        elif validators.url(mane_input):
            mane_summary_source = resolve_mane_summary_url(mane_input)
            logger.info("Resolved MANE summary URL: %s", mane_summary_source)
        else:
            sys.exit(f"\nERROR: MANE input is invalid: {mane_input}. Provide a valid local file path or URL")

        mane_file = resolve_local_or_url(mane_summary_source, file_path, label='MANE summary')
        mane_class_df, mane_meta = load_mane_transcript_map(mane_file)

        resource_files_dict["mane_status"] = "mane_status.tsv.gz"

        write_table(mane_class_df,
                    out_dir=out_path,
                    filename=resource_files_dict['mane_status'],
                    sep='\t',
                    index=False)

    
    # LOGGING AND CONFIG FILE WRITE

    logger.info("Resource bundle written to: %s", out_path)
    # Log the resource files create
    log_file = os.path.join(out_path, 'resource_metadata.txt')
    with open(log_file, 'w') as (wh):
        import datetime
        current_datetime = datetime.datetime.now()
        formatted_datetime = current_datetime.strftime('%Y-%m-%d %H:%M:%S')
        wh.write(f"Time of creation: {formatted_datetime}\n")
        wh.write(f"\nSource of data: {gtf_input}\n")
        wh.write(f"Genome build: {args.genome}\n")
        wh.write(f"\nAnnotation summary source: {annot_summary_input}\n")
        if mane_enabled:
            wh.write(f"MANE summary source (input): {mane_input}\n")
            wh.write(f"MANE summary source (resolved): {mane_summary_source}\n")
        else:
            wh.write("MANE annotation resource: disabled and not used\n")
        wh.write(f"bundle directory name: {Path(out_path).name}\n")

        wh.write("\nInput file checksums (SHA256):\n")
        for label, maybe_path in [
            ("gtf", gtf_file),
            ("annotation_summary", annot_summary_file),
            ("mane_summary", mane_file),
        ]:
            if maybe_path and os.path.isfile(maybe_path):
                wh.write(f"{label}: {os.path.basename(maybe_path)}\t{file_sha256(maybe_path)}\n")


        wh.write("\nGenerated resource files:\n")
        for k, v in resource_files_dict.items():
            target_path = os.path.join(out_path, v)
            if os.path.isfile(target_path):
                wh.write(f"{k}: {v}\tsha256={file_sha256(target_path)}\n")
            else:
                wh.write(f"{k}: {v}\n")
    
    write_config_file(out_path, resource_files_dict)


def drop_dup_col(df):
    cols_to_drop = []
    for i, col1 in enumerate(df.columns):
        for j, col2 in enumerate(df.columns):
            if i < j and df.iloc[:, i].equals(df.iloc[:, j]):
                cols_to_drop.append(j)

    cols_to_keep = [i for i in range(df.shape[1]) if i not in cols_to_drop]
    df = df.iloc[:, cols_to_keep]
    dup_cols = df.columns[df.columns.duplicated(keep=False)].unique()
    if dup_cols.any():
        for col_name in dup_cols:
            duplicate_col_positions = [i for i, name in enumerate(df.columns) if name == col_name]
            non_na_counts = [df.iloc[:, pos].count() for pos in duplicate_col_positions]
            pos_to_drop = duplicate_col_positions[non_na_counts.index(min(non_na_counts))]
            cols_to_keep = [i for i in range(df.shape[1]) if i != pos_to_drop]
            df = df.iloc[:, cols_to_keep]

    return df


def sort_bed_df(df):
    
    df = df.replace({'chr': ''}, regex=True)
    df = df.loc[pd.to_numeric((df['chrom']), errors='coerce').sort_values().index]
    df['chrom'] = 'chr' + df['chrom'].astype(str)
    df = df.groupby(['chrom'], sort=False, group_keys=True).apply((lambda x: x.sort_values(['start', 'end'], ascending=True))).reset_index(drop=True)
    return df


def write_config_file(out_path, resource_files_dict):
    config_file = os.path.join(out_path, 'manifest.yaml')
    with open(config_file, 'w') as wh:
        wh.write("resource_files:\n")
        for k, v in resource_files_dict.items():
            wh.write(f"  {k}: {v}\n")
        wh.write("\nrequired_keys:\n")
        for k in resource_files_dict.keys():
            wh.write(f"  - {k}\n")


def file_sha256(path: str | os.PathLike[str]) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_table(
    df: pd.DataFrame,
    *,
    out_dir: Path,
    filename: str,
    sep: str = "\t",
    index: bool = False,
):
    """
    Write a DataFrame to a compressed file whose name comes from the resource dict.

    Parameters
    ----------
    df : pandas.DataFrame
    out_dir : str or Path
        Directory where the file will be written.
    filename : str
        Filename from resource_files_dict (no path).
        Supports .zip and .gz.
    """
    out_dir = Path(out_dir)
    out_path = out_dir / filename
   
    df.to_csv(
            out_path,
            sep=sep,
            index=index,
            compression="gzip",
        )

if __name__ == '__main__':
    args = parse_args()
    main(args)
