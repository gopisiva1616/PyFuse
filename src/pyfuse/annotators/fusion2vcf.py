#!/usr/bin/env python3
"""
Convert a fusion DataFrame to a BND VCF file with alphabetic IDs and debugging support.
"""


import os
import argparse
import logging
import pandas as pd
import string


def get_vcf_header_template_path():
    """
    Returns the default vcfheader template path in the config directory.
    """
    this_dir = os.path.dirname(os.path.abspath(__file__))
    conf_path = os.path.abspath(os.path.join(this_dir, '../../config'))
    return os.path.join(conf_path, 'vcfheader.template.txt')

def get_bnd_ids(idx):
    """
    Returns (id1, id2) for a given fusion index: a/b, c/d, ... z/Z, aa/ab, ...
    Handles arbitrarily large numbers by using repeated letters: aaa, aab, ...
    """
    alphabet = list(string.ascii_lowercase) + list(string.ascii_uppercase)
    n = len(alphabet)
    def letter(num):
        # Converts a number to a base-n string using alphabet
        s = ''
        while True:
            s = alphabet[num % n] + s
            num = num // n
            if num == 0:
                break
        return s
    i1 = idx * 2
    i2 = idx * 2 + 1
    id1 = f"bnd_{letter(i1)}_bnd_{letter(i2)}"
    id2 = f"bnd_{letter(i2)}_bnd_{letter(i1)}"
    return id1, id2

def extract_transcript(annotation):
    # Returns transcript ID from annotation string
    if isinstance(annotation, str) and '|' in annotation:
        return annotation.split('|')[-1]
    return '.'

def fusion_df_to_bnd_vcf(df, vcf_path, vcf_header_path=None, debug_path=None, logger=None):
    # Use config-based template path if not provided
    if vcf_header_path is None:
        vcf_header_path = get_vcf_header_template_path()
    if vcf_header_path and os.path.exists(vcf_header_path):
        with open(vcf_header_path) as f:
            vcf_header = [line.rstrip() for line in f if line.strip() and not line.strip().startswith('## VCF Header Template')]
    else:
        vcf_header = [
            '##fileformat=VCFv4.5',
            '##source=PyFuse',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
        ]
    gene_col = "5'-3'Gene_Partners"
    debug_rows = []
    with open(vcf_path, 'w') as vcf:
        for line in vcf_header:
            vcf.write(line + '\n')
        for idx, row in df.iterrows():
            chrom5, pos5 = str(row["5'co-ordinate"]).split(':')
            chrom3, pos3 = str(row["3'co-ordinate"]).split(':')
            id5, id3 = get_bnd_ids(idx)
            alt5 = f"N]{chrom3}:{pos3}]"
            alt3 = f"N]{chrom5}:{pos5}]"
            transcript_5p = extract_transcript(row.get("5'_Exon_Annotation", "."))
            transcript_3p = extract_transcript(row.get("3'_Exon_Annotation", "."))
            transcript_field = f"{transcript_5p},{transcript_3p}" if transcript_5p != '.' and transcript_3p != '.' else transcript_5p if transcript_3p == '.' else transcript_3p
            fusion_id = row.get('Fusion_id', idx+1)
            info_common = [
                f"MATEID={id3}",
                "SVTYPE=BND",
                "EVTYPE=Fusion",
                f"GENE={row.get(gene_col, '.')}",
                f"FSN_LOC={row.get('Fusion_Annotation', '.')}",
                f"FRM_STAT={row.get('Frame_Status', '.')}",
                f"TRANSCRIPT={transcript_field}",
                f"FSN_ANNOT={row.get('Blacklist_Annotation', '.')}",
                f"MATEGENE={row.get('Fusion_Position', '.')}",
                f"HISTOLOGY={row.get('Histology', '.')}",
                f"PRESENT_IN_COSMIC={row.get('Present_in_COSMIC', '.')}",
                f"PRESENT_IN_GTEX={row.get('Present_in_GTEX', '.')}",
                f"AVERAGE_EXPRESSION={row.get('Average_Expression', '.')}",
                f"NUMBER_OF_TISSUES={row.get('Number_of_Tissues_that_contain_fusion', '.')}"
            ]
            info_common.append(f"FUSION_ID={fusion_id}")
            vcf.write(f"{chrom5}\t{pos5}\t{id5}\tN\t{alt5}\t.\t.\t" + ";".join(info_common) + "\n")
            info_common[0] = f"MATEID={id5}"
            vcf.write(f"{chrom3}\t{pos3}\t{id3}\tN\t{alt3}\t.\t.\t" + ";".join(info_common) + "\n")
            if debug_path:
                debug_rows.append({
                    'fusion_idx': idx,
                    'id5': id5,
                    'id3': id3,
                    'chrom5': chrom5,
                    'pos5': pos5,
                    'chrom3': chrom3,
                    'pos3': pos3,
                    'alt5': alt5,
                    'alt3': alt3,
                    'info5': ";".join(info_common),
                })
    if debug_path and debug_rows:
        pd.DataFrame(debug_rows).to_csv(debug_path, sep='\t', index=False)
        if logger:
            logger.info(f"Debug info written to {debug_path}")

def main():
    parser = argparse.ArgumentParser(description="Convert fusion table to BND-style VCF with debug support.")
    parser.add_argument('-i', '--input', required=True, help='Input fusion file (Excel or TSV)')
    parser.add_argument('-o', '--output', required=True, help='Output VCF file')
    parser.add_argument('-c', '--vcf_header', required=False, help='VCF header template file')
    parser.add_argument('-d', '--debug', action='store_true', help='Enable debug output')
    parser.add_argument('--debug_path', required=False, help='Debug TSV output path')
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.debug else logging.INFO)
    logger = logging.getLogger("fusion_df_to_bnd_vcf")

    # Read input
    if args.input.endswith('.xlsx'):
        df = pd.read_excel(args.input)
    else:
        df = pd.read_csv(args.input, sep=None, engine='python')

    debug_path = args.debug_path if args.debug else None
    fusion_df_to_bnd_vcf(df, args.output, vcf_header_path=args.vcf_header, debug_path=debug_path, logger=logger)
    logger.info(f"VCF written to {args.output}")

if __name__ == "__main__":
    main()
