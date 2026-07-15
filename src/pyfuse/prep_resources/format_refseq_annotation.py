#!/usr/bin/env python3

"""

Script to parse and annotate refseq_gene model downloaded from the UCSC
table browser.

https://genome.ucsc.edu/cgi-bin/hgTables


Date: 14/06/2021

"""

import sys
import pandas as pd

refseq_table_file = sys.argv[1]
output_file = sys.argv[2]


chrom_list = ['chr'+str(i) for i in list(range(1, 23)) + ['X', 'Y']]
req_cols = ['chrom', 'strand', 'exonStarts', 'exonEnds', 'exonCount',
            'txStart', 'txEnd', 'name2', 'name']
raw_refseq_df = pd.read_csv(refseq_table_file, sep='\t', usecols=req_cols)
raw_refseq_df.rename(columns={'name': 'transcript',
                              'name2': 'gene'}, inplace=True)

# removing genes from that are not in 24chroms
refseq_df = raw_refseq_df[raw_refseq_df['chrom'].isin(chrom_list)]

# selecting longest transcript for each gene

refseq_df = refseq_df.assign(transcript_length=refseq_df['txEnd'] - refseq_df['txStart'])
#refseq_df.apply(lambda x: diff(x['txEnd'], x['txStart']), axis=1)
refseq_df = refseq_df.loc[refseq_df.groupby(['gene'], sort=False)['transcript_length'].idxmax()].reset_index(drop=True)

#split by exons and adding exon numbers
refseq_exploded = refseq_df[['exonStarts', 'exonEnds']].apply(lambda x: [v.split(',') for v in x]).apply(pd.Series.explode)
refseq_df = pd.concat([refseq_df[['transcript', 'chrom', 'strand',  'txStart',
                                  'txEnd', 'gene']], refseq_exploded], axis=1)
refseq_df['exon_number'] = refseq_df.groupby(['transcript', 'chrom', 'strand',
                                              'txStart', 'txEnd', 'gene']).cumcount()+1
refseq_df = refseq_df[['chrom', 'txStart', 'txEnd', 'strand', 'exonStarts', 'exonEnds', 'exon_number', 'gene', 'transcript']]

refseq_df.to_csv(output_file, sep='\t', index=False, header=False)
