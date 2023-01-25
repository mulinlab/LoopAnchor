#!/usr/bin/env python
r"""
This script generates training dataset for DeepAnchor.
Please include following data within a working directory and arrange them with structure as shown below:
work_dir
    ----raw
        ----target.bed           # Targeted chromosome intervals in bed format.
        ----CTCF_peak.bed.gz     # The ChIP-seq peak files of CTCF.
        ----cbs.tsv       # Position of all CTCF binding sites. 
        ----dna_feature.npz      # One-hot representation of DNA sequence for CTCF binding sites.
        ----cadd_feature.npz     # Cadd features of DNA sequence for CTCF binding sites.

You should prepare target.bed, CTCF-peak.bed.gz for specific sample.

cbs.tsv, dna_feature.npz, cadd_feature.npz are non-cell-type-specific and can be downloaded from tutorial.

Usage
-----
python DeepAnchor_input.py work_dir

"""

import os
import sys
import pandas as pd
import numpy as np
import tensorflow as tf


bedpe_columns = ['chrom1','start1','end1','chrom2','start2','end2','name','score','strand1','strand2']
peak_columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'summit']
cbs_columns = ['chrom','start','end','strand','score']


def generate_marked_cbs(file_cbs, file_target):
    r"""
    Use chromosome intervals to mark CTCF binding sites.

    1. file_target should not contain overlaped intervals.
    2. We mark file_cbs with file_target and generate three new columns for file_cbs:
                    in_bed: 1 if cbs is contained in any intervals in file_target else 0
        exclusively_in_bed: 1 if cbs is exclusively contained in any intervals in file_target else 0
            maximum_in_bed: 1 if cbs is the one with highest score in any intervals in file_target esle 0
            
    Parameters
    ----------
    file_cbs
        File name of CTCF binding sites.
    file_targed
        File name of targeted chromosome intervals.


    Return
    ------
    df_cbs
        CBS data frame with 3 new columns.
    """
    
    df_cbs = pd.read_csv(file_cbs, sep='\t', names=cbs_columns)

    os.system('bedtools intersect -a {} -b {} -loj > temp.intersect'.format(file_cbs, file_target))
    bed_columns = ['bed_chrom','bed_start','bed_end']
    df_matrix = pd.read_csv('temp.intersect', sep='\t', names=cbs_columns + bed_columns)
    
    df_matrix = df_matrix.sort_values(['score'], ascending=False)
    
    list1, list2, list3 = [], [], []
    for bed, sub_cbs in df_matrix.groupby(['bed_chrom','bed_start','bed_end']):
        if bed[0] == '.':
            pass
        elif sub_cbs.shape[0] > 1:
            list1 += sub_cbs.index.tolist()
            list3.append(sub_cbs.index.tolist()[0])
        else:
            list1 += sub_cbs.index.tolist()
            list2 += sub_cbs.index.tolist()
            list3 += sub_cbs.index.tolist()

    df_cbs['in_bed'] = 0
    df_cbs['exclusively_in_bed'] = 0
    df_cbs['maximum_in_bed'] = 0
    df_cbs.loc[list1, 'in_bed'] = 1
    df_cbs.loc[list2, 'exclusively_in_bed'] = 1
    df_cbs.loc[list3, 'maximum_in_bed'] = 1

    return df_cbs


def main(argv=sys.argv):
    work_dir = sys.argv[1]

    os.chdir(os.path.expanduser(work_dir))

    if not os.path.exists('DeepAnchor'):
        os.mkdir('DeepAnchor')
    
    file_target = './raw/target.bed'
    file_peak = './raw/CTCF_peak.bed.gz'
    file_cbs = './raw/cbs.tsv'
    file_dna_features = './raw/dna_feature.npz'
    file_cadd_features = './raw/cadd_feature.npz'
        
    print('>>> Preprocess ChIP-seq peaks')
    if not os.path.exists('./DeepAnchor/CTCF_peaks.bed'):
        df_peak = pd.read_csv(file_peak, sep='\t', names=peak_columns)
        df_peak = df_peak.sort_values(['chrom','start'])
        df_peak.to_csv('temp.peak.bed', sep='\t', index=False, header=False)
        os.system('bedtools merge -i temp.peak.bed > ./DeepAnchor/CTCF_peaks.bed')


    print('>>> Mark cbs with targets and CTCF ChIP-seq peaks')
    if not os.path.exists('./DeepAnchor/marked_cbs.tsv'):
        df_anchor_mark_cbs = generate_marked_cbs(file_cbs, file_target)
        df_peak_mark_cbs = generate_marked_cbs(file_cbs, './DeepAnchor/CTCF_peaks.bed')
        df_marked_cbs = df_anchor_mark_cbs.merge(df_peak_mark_cbs, on=['chrom','start','end','strand','score'])
        df_marked_cbs.columns = ['chrom','start','end','strand','score','in_anchor','exclusively_in_anchor','maximum_in_anchor','in_peak','exclusively_in_peak','maximum_in_peak']
        df_marked_cbs.to_csv('./DeepAnchor/marked_cbs.tsv', sep='\t', index=False)


    print('>>> Select positive and negative data')
    df_marked_cbs = pd.read_csv('./DeepAnchor/marked_cbs.tsv', sep='\t')
    df_marked_cbs['selected'] = 0
    df_pos = df_marked_cbs[(df_marked_cbs['exclusively_in_anchor'] == 1) & 
                             (df_marked_cbs['in_peak'] == 1) ]
    pos_index = df_pos.index.tolist()

    # df_neg = df_marked_cbs[(df_marked_cbs['in_anchor'] ==0) &
    #                          (df_marked_cbs['in_peak'] == 1)]
    df_neg = df_marked_cbs[df_marked_cbs['in_anchor'] == 0]

    selected_neg_index = list(np.random.choice(df_neg.index.tolist(), len(pos_index)))
    df_marked_cbs.loc[pos_index+selected_neg_index, 'selected'] = 1


    print('>>> Add features')
    if not os.path.exists('./DeepAnchor/train.npz') \
       or not os.path.exists('./DeepAnchor/test.npz') \
       or not os.path.exists('./DeepAnchor/valid.npz') :

        print('> load dna and cadd features')
        dna=np.load(file_dna_features)['dna'].astype(np.float32)
        cadd_features=np.load(file_cadd_features)['cadd_features'].astype(np.float32)
    
    
        print('> select training dataset')
        df_total = df_marked_cbs[df_marked_cbs['selected'] == 1]
        df_total['label'] = 0
        df_total.loc[pos_index, 'label'] = 1
        df_train = df_total[~df_total['chrom'].isin(['chr17','chr18','chr19','chr20','chr21','chr22','chrX'])]
        df_valid = df_total[df_total['chrom'].isin(['chr17','chr18'])]
        df_test = df_total[df_total['chrom'].isin(['chr19','chr20','chr21','chr22','chrX'])]
    
        print('Number of train, valid, test data:')
        print('train(chr1-16): {}'.format(df_train.shape[0]))
        print('valid(chr17,chr18): {}'.format(df_valid.shape[0]))
        print('test(chr19-X): {}'.format(df_test.shape[0]))
    
        train_index = df_train.index.tolist()
        train_dna = dna[train_index]
        train_cadd = cadd_features[train_index]
        train_Features = np.concatenate([train_cadd, train_dna], axis=2)
        train_label = df_total.loc[train_index, 'label'].values
            
        valid_index = df_valid.index.tolist()
        valid_dna = dna[valid_index]
        valid_cadd = cadd_features[valid_index]
        valid_Features = np.concatenate([valid_cadd, valid_dna], axis=2)
        valid_label = df_total.loc[valid_index, 'label'].values

        test_index = df_test.index.tolist()
        test_dna = dna[test_index]
        test_cadd = cadd_features[test_index]
        test_Features = np.concatenate([test_cadd, test_dna], axis=2)
        test_label = df_total.loc[test_index, 'label'].values

        print('> save data in npz')
        np.savez('./DeepAnchor/train.npz', Features=train_Features, label=train_label)
        np.savez('./DeepAnchor/valid.npz', Features=valid_Features, label=valid_label)
        np.savez('./DeepAnchor/test.npz', Features=test_Features, label=test_label)

    print('>>> generate total npz')
    if not os.path.exists('./DeepAnchor/total.npz'):
        dna=np.load(file_dna_features)['dna'].astype(np.float32)
        cadd_features=np.load(file_cadd_features)['cadd_features'].astype(np.float32)

        total_Features = np.concatenate([cadd_features, dna], axis=2)
        np.savez('./DeepAnchor/total.npz', Features=total_Features)
    
    os.system('rm temp.*')

if __name__ == '__main__':
    main()
