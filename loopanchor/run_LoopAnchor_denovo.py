#!/usr/bin/env python
r"""
Use Loop extrusion model to do loop prediction.
Please include following data within a work_dir and arrange them like that:
work_dir
    ----raw
        ----CTCF_peak.bed.gz     # The ChIP-seq peak files of CTCF
run: run_LoopAnchor_denovo.py  work_dir
"""

import sys 
import os

import pandas as pd
import numpy as np

from tqdm import tqdm
from src.model import LoopAnchor

from pathlib import Path

work_dir = sys.argv[1]


dir_input = Path(work_dir, 'raw')
dir_output = Path(work_dir, 'LoopAnchor')

os.makedirs(dir_output, exist_ok=True)

# load peak file
# input files
f_peak = dir_input / 'CTCF_peak.bed.gz'
f_cbs = dir_input / 'scored_cbs.tsv'

peak_columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'summit']
bedpe_columns = ['chrom1','start1','end1','chrom2','start2','end2','name', 'score','strand1','strand2']
df_peak = pd.read_csv(f_peak, sep='\t', names=peak_columns)
df_peak = df_peak.sort_values(['chrom','start','end'])


cbs_columns = ['chrom','start','end','strand','score', 'anchor_score']
df_cbs = pd.read_csv(f_cbs, sep='\t', names=cbs_columns)
df_cbs = df_cbs.sort_values('anchor_score', ascending=False)



print('>>> Add cbs on peak')
list_peak_cbs = []

for chrom in tqdm(df_peak['chrom'].unique()):
    sub_peak = df_peak[df_peak['chrom'] == chrom]
    sub_cbs = df_cbs[df_cbs['chrom'] == chrom]
    for _, row in sub_peak.iterrows():
        chrom, start, end, h = row[['chrom', 'start','end','signalValue']]
        select_cbs = sub_cbs[(sub_cbs['start'] > start) & (sub_cbs['end'] < end)]
        if len(select_cbs) == 0:
            list_peak_cbs.append([chrom, start, end, h, -1, -1, -1])
        else:
            strand, score, anchor = select_cbs.iloc[0, 3:]
            list_peak_cbs.append([chrom, start, end, h, strand, score, anchor])

df_peak_cbs = pd.DataFrame(list_peak_cbs, columns=['chrom','start','end','h','strand', 'score', 'anchor'])
df_peak_cbs = df_peak_cbs[df_peak_cbs['score'] != -1]
df_peak_cbs['index'] = np.arange(len(df_peak_cbs)) + 1


print('>>> Generate all loops')
list_total = []
list_loop = []
for chrom in tqdm(df_peak_cbs['chrom'].unique()):
    sub_df1 = df_peak_cbs[df_peak_cbs['chrom'] == chrom]
    for _, row1 in sub_df1.iterrows():
        chrom1, start1, end1, h1, strand1, score1, anchor1, index1 = row1
        sub_df2 = sub_df1[(sub_df1['start'] - start1 > 10000) & (sub_df1['start'] - start1 < 1000000)]
        for index2, row2 in sub_df2.iterrows():
            chrom2, start2, end2, h2, strand2, score2, anchor2, index2 = row2
            distance = start2 - start1
            label = 0 
            count = 0
            list_total.append([index1,index2,h1,h2,distance,strand1,score1,anchor1,strand2,score2,anchor2,label,count])
            list_loop.append([chrom1, start1, end1, chrom2, start2, end2, strand1, strand2])

df_total = pd.DataFrame(list_total, columns = ['index1' ,'index2', 'h1', 'h2', 'distance', 'strand1', 'score1','anchor1','strand2','score2','anchor2','label','count'])
df_loop = pd.DataFrame(list_loop, columns = ['chrom1','start1','end1','chrom2','start2','end2','strand1','strand2'])
df_loop['name'] = '.'
df_loop['score'] = 0


print('>>> Make prediction for all loops')
df_LoopAnchor = LoopAnchor(df_peak_cbs, df_total[['index1', 'index2', 'h1', 'h2', 'distance', 'strand1', 'anchor1', 'strand2', 'anchor2']])
df_loop['LoopAnchor'] = df_LoopAnchor['pred']

df_output = df_loop[bedpe_columns+['LoopAnchor']]
df_output.to_csv(dir_output / 'LoopAnchor_pred.bedpe', sep='\t', index=False, header=False)
