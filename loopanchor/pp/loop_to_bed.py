#!/usr/bin/env python
r"""
Extract anchor (bed format) from loop file (bedpe format)
"""


import os
import pandas as pd

import sys

f_in = sys.argv[1]
f_out = sys.argv[2]

bedpe_columns = ['chrom1','start1','end1','chrom2','start2','end2','name','score','strand1','strand2']

df_loop = pd.read_csv(f_in, sep='\t', names=bedpe_columns)
anchor1=pd.DataFrame(df_loop[['chrom1','start1','end1']].values, columns=['chrom','start','end'])
anchor2=pd.DataFrame(df_loop[['chrom2','start2','end2']].values, columns=['chrom','start','end'])
total=pd.concat([anchor1, anchor2]).sort_values(['chrom','start','end'])
total.to_csv('temp.total.anchor', sep='\t', index=False, header=False)
os.system('bedtools merge -i temp.total.anchor > {}'.format(f_out))

os.system('rm temp.total.anchor')
