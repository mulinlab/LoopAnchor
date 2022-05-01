from tqdm import tqdm
import pandas as pd

def bed_in_pos(df_bed, pos, uniq_chrom=False):
    chrom, start, end = pos

    if uniq_chrom == True:
        sub_bed = df_bed
    else:
        sub_bed = df_bed[df_bed['chrom'] == chrom]

    sub_bed = sub_bed[(sub_bed['start'] >= start) & (sub_bed['end'] <= end)]

    if len(sub_bed) == 0:
        index = -1

    else:
        index = sub_bed['index'].iat[0]

    return index

def motif_on_chip(df_peak, df_motif):
    list_motif_index = []
    for chrom in tqdm(df_peak['chrom'].unique()):
        sub_peak = df_peak[df_peak['chrom'] == chrom]
        sub_motif = df_motif[df_motif['chrom'] == chrom]
        for _, row in sub_peak.iterrows():
            chrom, start, end = row[['chrom','start','end']]
            motif_index = bed_in_pos(sub_motif, [chrom, start, end], uniq_chrom=True)
            list_motif_index.append(motif_index)
    return list_motif_index


def pos_in_bed(df_bed, pos, uniq_chrom=False):
    chrom, start, end = pos

    if uniq_chrom == True:
        sub_bed = df_bed
    else:
        sub_bed = df_bed[df_bed['chrom'] == chrom]

    sub_bed = sub_bed[(sub_bed['start'] < start) & (sub_bed['end'] > end)]

    if len(sub_bed) == 0:
        index = -1

    else:
        index = sub_bed['index'].iat[0]

    return index



def peak_on_loop(df_loop, df_peak):
    dict_peak = {chrom: df_peak[df_peak['chrom'] == chrom] for chrom in df_peak['chrom'].unique()}

    left = df_loop[['chrom1','start1','end1']]
    left_peak_index = []

    for _, row in tqdm(left.iterrows()):
        chrom, start, end = row
        if chrom not in dict_peak:
            peak_index = -1
        else:
            peak_index =bed_in_pos(dict_peak[chrom], [chrom, start, end], uniq_chrom=True)
        left_peak_index.append(peak_index)

    right = df_loop[['chrom2','start2','end2']]
    right_peak_index = []

    for _, row in tqdm(right.iterrows()):
        chrom, start, end = row
        if chrom not in dict_peak:
            peak_index = -1
        else:
            peak_index = bed_in_pos(dict_peak[chrom], [chrom, start, end], uniq_chrom=True)
        right_peak_index.append(peak_index)
    
    return left_peak_index, right_peak_index
