import pandas as pd
from math import exp

def LoopAnchor(df_peak_cbs, df_loop):
    r"""
    Model to predict Cohesin-mediated loops.

    Parameters
    ----------

    df_peak_cbs
        CTCF ChIP-seq peaks (marked with cbs).
    df_loop
        All positive loops mediated by Cohesin/CTCF.
    """

    # constants
    w_ori = 3
    w_ctcf = 8.5
    lambda_ = 3000000
    
    
    #### load peak file
    n = len(df_peak_cbs)
    ctcf = df_peak_cbs[['h', 'anchor']].apply(lambda s: s[0] * s[1], axis=1).mean()

    
    peak_dict = {}
    for index, row in df_peak_cbs.iterrows():
        chrom,start,end,h,strand,cbs_score,score,index = row
        peak_dict[index] = [chrom, (start+end)/2, h, strand, cbs_score, score]
    
    #### load training file
    y_true = []
    y_pred = []
    
    list_pred = []
    list_processicivity = []
    list_lc = []
    list_pij = []
    list_ori = []
    
    for index, row in df_loop.iterrows():
        index1, index2, h1, h2, dist, strand1, score1, strand2, score2 = row
        strand1 = -1 if strand1 == '-' else 1
        strand2 = -1 if strand2 == '-' else 1
        
        h1 = h1 * score1
        h2 = h2 * score2

        lc = 1
    
        p_ij = (h1 / (h1 + ctcf * w_ctcf)) * (h2/(h2 + ctcf * w_ctcf))
        ori = w_ori ** ((strand1 - strand2 -2)/2)
        processicivity = exp(-dist/lambda_)
    
    
        for k in range(index1+1, index2):
            lc = lc * (1 - peak_dict[k][2]/(peak_dict[k][2] + ctcf * w_ctcf))
    
        prediction = ori * p_ij * lc
    
        list_pij.append(p_ij)
        list_ori.append(ori)
        list_pred.append(prediction)
        list_lc.append(lc)
        list_processicivity.append(processicivity)
    
    df_training = df_loop.copy(deep=True)
    df_training['p_ij'] = list_pij
    df_training['ori'] = list_ori
    df_training['processicivity'] = list_processicivity
    df_training['pred'] = list_pred
    df_training['lc'] = list_lc

    return df_training
