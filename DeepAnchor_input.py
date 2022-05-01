#!/usr/bin/env python
'''
This script generates training dataset for DeepAnchor.
Please include following data within a work_dir and arrange them like that:
work_dir
    ----raw
        ----loop.bedpe           # ChIA-PET or other types of loop files in bedpe format
        ----CTCF_peak.bed.gz     # The ChIP-seq peak files of CTCF
        ----CTCF_motif.tsv       # position of all CTCF binding sites 
        ----dna_feature.npz      # one-hot representation of DNA sequence for CTCF binding sites
        ----cadd_feature.npz     # cadd features of DNA sequence for CTCF binding sites

You should prepare loop.bedpe, CTCF-peak.bed.gz for specific sample, while CTCF_motif.tsv, 
dna_feature.npz, cadd_feature.npz are non-cell-type-specific and can be downloaded from tutorial.

usage: python DeepAnchor_input.py work_dir
'''

import os
import sys
import pandas as pd
import numpy as np
import tensorflow as tf


bedpe_columns = ['chrom1','start1','end1','chrom2','start2','end2','name','score','strand1','strand2']
peak_columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'summit']
motif_columns = ['chrom','start','end','strand','score']


def generate_bed_mark_motif(file_motif, file_bed):
    '''
    1. file_bed should contain no overlaps
    2. We mark file_motif with file_bed and generate three marks:
                    in_bed: 1 if motif is contained in any intervals in file_bed else 0
        exclusively_in_bed: 1 if motif is exclusively contained in any intervals in file_bed else 0
            maximum_in_bed: 1 if motif is the one with highest score in any intervals in file_bed esle 0
    '''
    
    df_motif = pd.read_csv(file_motif, sep='\t', names=motif_columns)

    os.system('bedtools intersect -a {} -b {} -loj > temp.motif'.format(file_motif, file_bed))
    bed_columns = ['bed_chrom','bed_start','bed_end']
    df_matrix = pd.read_csv('temp.motif', sep='\t', names=motif_columns + bed_columns)
    
    df_matrix = df_matrix.sort_values(['score'], ascending=False)
    
    list1, list2, list3 = [], [], []
    for bed, sub_motif in df_matrix.groupby(['bed_chrom','bed_start','bed_end']):
        if bed[0] == '.':
            pass
        elif sub_motif.shape[0] > 1:
            list1 += sub_motif.index.tolist()
            list3.append(sub_motif.index.tolist()[0])
        else:
            list1 += sub_motif.index.tolist()
            list2 += sub_motif.index.tolist()
            list3 += sub_motif.index.tolist()

    df_motif['in_bed'] = 0
    df_motif['exclusively_in_bed'] = 0
    df_motif['maximum_in_bed'] = 0
    df_motif.loc[list1, 'in_bed'] = 1
    df_motif.loc[list2, 'exclusively_in_bed'] = 1
    df_motif.loc[list3, 'maximum_in_bed'] = 1

    return df_motif


def main(argv=sys.argv):
    work_dir = sys.argv[1]

    os.chdir(os.path.expanduser(work_dir))

    if not os.path.exists('DeepAnchor'):
        os.mkdir('DeepAnchor')
    
    file_bedpe = './raw/loop.bedpe'
    file_peak = './raw/CTCF_peak.bed.gz'
    file_motif = './raw/CTCF_motif.tsv'
    file_dna_onehot = './raw/dna_feature.npz'
    file_cadd_features = './raw/cadd_feature.npz'


    print('>>> generate total anchors')
    if not os.path.exists('./DeepAnchor/total_anchors.bed'):
        df_loop = pd.read_csv(file_bedpe, sep='\t', names=bedpe_columns)
        anchor1=pd.DataFrame(df_loop[['chrom1','start1','end1']].values, columns=['chrom','start','end'])
        anchor2=pd.DataFrame(df_loop[['chrom2','start2','end2']].values, columns=['chrom','start','end'])
    
        total=pd.concat([anchor1, anchor2]).sort_values(['chrom','start','end'])
        total.to_csv('temp.total.anchor', sep='\t', index=False, header=False)
        os.system('bedtools merge -i temp.total.anchor > ./DeepAnchor/total_anchors.bed')

        
    print('>>> generate file peaks')
    if not os.path.exists('./DeepAnchor/CTCF_peaks.bed'):
        df_peak = pd.read_csv(file_peak, sep='\t', names=peak_columns)
        df_peak = df_peak.sort_values(['chrom','start'])
        df_peak.to_csv('temp.peak.bed', sep='\t', index=False, header=False)
        os.system('bedtools merge -i temp.peak.bed > ./DeepAnchor/CTCF_peaks.bed')


    print('>>> mark motif with anchor and peak')
    if not os.path.exists('./DeepAnchor/marked_motif.tsv'):
        df_anchor_mark_motif = generate_bed_mark_motif(file_motif, './DeepAnchor/total_anchors.bed')
        df_peak_mark_motif = generate_bed_mark_motif(file_motif, './DeepAnchor/CTCF_peaks.bed')
        df_marked_motif = df_anchor_mark_motif.merge(df_peak_mark_motif, on=['chrom','start','end','strand','score'])
        df_marked_motif.columns = ['chrom','start','end','strand','score','in_anchor','exclusively_in_anchor','maximum_in_anchor','in_peak','exclusively_in_peak','maximum_in_peak']
        df_marked_motif.to_csv('./DeepAnchor/marked_motif.tsv', sep='\t', index=False)


    print('>>> Select positive and negative data')
    df_marked_motif = pd.read_csv('./DeepAnchor/marked_motif.tsv', sep='\t')
    df_marked_motif['selected'] = 0
    df_pos = df_marked_motif[(df_marked_motif['exclusively_in_anchor'] == 1) & 
                             (df_marked_motif['in_peak'] == 1) ]
    pos_index = df_pos.index.tolist()

    # df_neg = df_marked_motif[(df_marked_motif['in_anchor'] ==0) &
    #                          (df_marked_motif['in_peak'] == 1)]
    df_neg = df_marked_motif[df_marked_motif['in_anchor'] == 0]

    selected_neg_index = list(np.random.choice(df_neg.index.tolist(), len(pos_index)))
    df_marked_motif.loc[pos_index+selected_neg_index, 'selected'] = 1


    print('>>> Add features')
    if not os.path.exists('./DeepAnchor/train.npz') \
       or not os.path.exists('./DeepAnchor/test.npz') \
       or not os.path.exists('./DeepAnchor/valid.npz') :

        print('> load dna and cadd features')
        dna=np.load(file_dna_onehot)['dna'].astype(np.float32)
        cadd_features=np.load(file_cadd_features)['cadd_features'].astype(np.float32)
    
    
        print('> select training dataset')
        df_total = df_marked_motif[df_marked_motif['selected'] == 1]
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
        dna=np.load(file_dna_onehot)['dna'].astype(np.float32)
        cadd_features=np.load(file_cadd_features)['cadd_features'].astype(np.float32)

        total_Features = np.concatenate([cadd_features, dna], axis=2)
        np.savez('./DeepAnchor/total.npz', Features=total_Features)
    
    os.system('rm temp.*')

if __name__ == '__main__':
    main()
