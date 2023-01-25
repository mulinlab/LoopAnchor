#!/usr/bin/env python
r"""
DeepAnchor was a Deep Learning method to identify CTCF binding sites located at specific targets.

It has two different mode: training mode and prediction mode.

When mode is train, it trains a classifier with Conv2D.

When mode is predict, it uses the classifier to predict the score of all cbs.


Usage
-----
python DeepAnchor.py  work_dir mode
"""

import os
import sys

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.preprocessing import OneHotEncoder

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import models, layers


def ROC_AUC(y_raw, y_pred):
    r"""
    Plot ROC curve.

    Parameters
    ----------
    y_raw
        Ground truth of label.
    y_pred
        Prediction score of label.

    """
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(9,4), sharex=True, sharey=True)

    FPR, TPR, _ = metrics.roc_curve(y_raw, y_pred)
    auc = metrics.roc_auc_score(y_raw, y_pred)
    ax1.plot(FPR, TPR, lw=2, label='ROC AUC = %0.2f' % auc)
    precision, recall, _ = metrics.precision_recall_curve(y_raw, y_pred)
    aupr = metrics.average_precision_score(y_raw, y_pred)
    ax2.plot(recall, precision, lw=2, label='AUPR = %0.2f' % aupr)
    ax1.legend()
    ax2.legend()
    ax1.set_xlabel('1-Specificity')
    ax1.set_ylabel('Sensitivity')
    ax2.set_xlabel('Recall')
    ax2.set_ylabel('Precision')
    ax1.set_box_aspect(1)
    ax2.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    

def plot_history(histories, key='binary_crossentropy'):
    r"""
    Plot the training history of DL model.
    """
    plt.figure(figsize=(8,6))
  
    for name, history in histories:
        val = plt.plot(history.epoch, history.history['val_'+key],
                     '--', label=name.title()+' Val')
        plt.plot(history.epoch, history.history[key], color=val[0].get_color(),
               label=name.title()+' Train')
  
    plt.xlabel('Epochs')
    plt.ylabel(key.replace('_',' ').title())
    plt.legend()
  
    plt.xlim([0,max(history.epoch)])
    plt.show()


def DeepAnchor(file_model, data):
    r"""
    
    Structure of DeepAnchor and training process.

    Parameters
    ----------
    file_model
        File to save the model.
    data
        Training data.

    return
        model
    """

    ### construct model
    model = models.Sequential()
    model.add(layers.Conv1D(32, 10, activation='relu', input_shape=(1000, 48)))
    model.add(layers.MaxPooling1D(10))
    model.add(layers.Dropout(0.5))
    model.add(layers.Conv1D(64, 10, activation='relu'))
    model.add(layers.MaxPooling1D(2))
    model.add(layers.Dropout(0.5))
    model.add(layers.Flatten())
    model.add(layers.Dense(512, activation='relu'))
    model.add(layers.Dense(2, activation='sigmoid'))
    model.summary()

    model.compile(optimizer='adam',
        loss='binary_crossentropy',
        metrics=['accuracy','binary_crossentropy', tf.keras.metrics.AUC()])
    
    earlystop_callback = keras.callbacks.EarlyStopping(
        monitor='val_binary_crossentropy', min_delta=0.0001,
        patience=3)

    ### loading data
    train_X, train_y = data['train_Features'], data['train_Label']
    valid_X, valid_y = data['valid_Features'], data['valid_Label']
    test_X, test_y = data['test_Features'], data['test_Label']
    
    ### training
    history_model = model.fit(train_X, 
        train_y, 
        epochs=20, 
        batch_size=50,
        validation_data=(valid_X, valid_y),
        callbacks=[earlystop_callback]
        )
    plot_history([('DeepAnchor', history_model)])
    model.save(file_model)

    ### testing
    pred_y = model.predict(test_X)[:, 1]
    ROC_AUC(test_y, pred_y)

    return model



def main(argv=sys.argv):
    work_dir = argv[1]
    mode = argv[2]


    # input and output
    os.chdir(os.path.expanduser(work_dir))
    file_cbs = './raw/cbs.tsv'
    file_scored_cbs = './scored_cbs.tsv'
    file_model = './DeepAnchor.model'

    # train
    enc = OneHotEncoder(categories='auto')
    if mode == 'train':
        ### load data
        train_data = np.load('./DeepAnchor/train.npz')
        train_X, train_y = train_data['Features'], train_data['label']
        train_y = enc.fit(train_y.reshape(-1,1)).transform(train_y.reshape(-1,1)).toarray()

        valid_data = np.load('./DeepAnchor/valid.npz')
        valid_X, valid_y = valid_data['Features'], valid_data['label']
        valid_y = enc.fit(valid_y.reshape(-1,1)).transform(valid_y.reshape(-1,1)).toarray()
        
        test_data = np.load('./DeepAnchor/test.npz')
        test_X, test_y = test_data['Features'], test_data['label']

        ### training
        data = {'train_Features': train_X, 'train_Label': train_y,
                'valid_Features': valid_X, 'valid_Label': valid_y,
                'test_Features': test_X, 'test_Label': test_y}
        model = DeepAnchor(file_model, data)
        
        pred_y = model.predict(test_X)[:, 1]
        df_test = pd.DataFrame({'test_label': test_y,
                                'LoopAnchor': pred_y})
                                
        df_test.to_csv('./DeepAnchor/test_result.tsv', sep='\t', index=False)


    elif mode == 'predict':
        df_cbs = pd.read_csv(file_cbs, sep='\t', names=['chrom','start','end','strand','score'])
        predict_X = np.load('./DeepAnchor/total.npz')['Features']

        model = keras.models.load_model(file_model)
        pred_y = model.predict(predict_X)[:,1]
        df_cbs['anchor_score'] = pred_y
        df_cbs.to_csv(file_scored_cbs, sep='\t', index=False, header=False)

    else:
        print('{} mode is not currently supported'.format(mode))


if __name__ == '__main__':
    main()
