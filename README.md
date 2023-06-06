# LoopAnchor

CCCTC-binding factor (CTCF) is a transcription regulator which is involved in many cellular processes. How CTCF recognizes DNA sequence to exert chromosome barrier or enhancer blocking effects remains to be fully interrogated. Despite many computational tools were developed to predict CTCF-mediated loops qualitatively or quantitatively, few could specially evaluate the regulatory potential of DNA sequence at CTCF binding sites (CBSs) and how it affects chromatin loop formation. Here, we developed a deep learning model, DeepAnchor, to precisely characterize the binding patterns for different types of CBSs. By incorporating base-wise genomic/epigenomic features, we revealed distinct chromatin and sequence features for CTCF-mediated insulation and looping at a high resolution, such as two sequence motifs flanking the core CTCF motif at loop-associated CBSs. Besides, we leveraged the predicted anchor score to optimize the loop extrusion model and achieved the best performance in predicting CTCF-anchored loops. We established a compendium of context-specific CTCF-anchored loops across 52 human tissue/cell types and found that genomic disruption of CTCF-anchored loops may represent a general causal mechanism of disease pathogenesis. These computational models, together with the established resource, could facilitate the mechanistic research on how the CTCF-mediated cis-regulatory elements (CREs) shapes context-specific gene regulation in cell development and disease progression.



<p align="center">
   <img src="https://github.com/mulinlab/LoopAnchor/blob/master/docs/source/flowchart.PNG?raw=True">
</p>

## Prepare data
To implement DeepAnchor, one should prepare four types of data:
1. Obtain CBSs by motif scan.
2. Download base-wise genomic/epigenomic features from [CADD v1.3](https://cadd.gs.washington.edu/download). Here we also provide preprocessed CADD features for CTCF binding sites.
3. CTCF ChIP-seq data for specific cell type (for example GM12878).
4. Cohesin ChIA-PET data for the same cell type. 


## Generate training data
### 1. Extract feature values
We totally use 44 genomic/epigenomic features from CADD database and 4 DNA sequence feature (A/T/G/C). For each CBS and each feature, ±500bp feature values are extracted. Therefore, the size of feature matrix for each CBS is 48 x 1000. Feature matrix is strandardized to facilitate downstream analyses. However, this step is very time consuming. So we provide the preprocessed feature matrix data which can be downloaded from here: [cadd feature matrix](http://www.mulinlab.org/LoopAnchor/cadd_feature.npz) and [dna feature matrix](http://www.mulinlab.org/LoopAnchor/dna_feature.npz).

### 2. Generate positive/negative CBSs
We need create a work_dir with structure as below and put five data files within it.
```
work_dir
    └── raw                   
        ├── CTCF_peak.bed.gz              
        ├── target.bed              
        ├── cbs.tsv
        ├── dna_feature.npz      
        └── cadd_feature.npz
```
More explainations of data requirement:

CTCF_peak.bed.gz: 
    (1) the ChIP-seq narrowPeak file for CTCF.
    (2) Columns: chrom,chromStart,chromEnd,name,score,strand,signalValue,pValue,qValue,summit.

target.bed:
    (1) Targeted chromosome intervals for CTCF binding sites.
    (2) Columns: chrom,start,end.

cbs.tsv:
    (1) position of all CTCF binding sites.
    (2) non-cell-type specific.
    (3) can be found in data folder.

dna_feature.npz:
    (1) one-hot representation of DNA sequence for CTCF binding sites.
    (2) non-cell-type specific.
    (3) downloaded as addressed in last section.

cadd_feature.npz:
    (1) cadd features of DNA sequence for CTCF binding sites
    (2) non-cell-type specific
    (3) downloaded as addressed in last section.

To generate P/N dataset, you can simply run following command:
```properties
python DeepAnchor_input.py work_dir
```


Tips:

If you have loop bedpe file, you can use ./pp/loop_to_bed.py to convert the loops to target.bed file.

```
python ./pp/loop_to_bed.py file_loop target_bed
```


DeepAnchr_input.py will generate an output folder named DeepAnchor in work_dir:
```
work_dir
    └── DeepAnchor  
        ├── CTCF_peaks.bed           # intersect CTCF_peak.bed.gz and cbs.tsv
        ├── marked_cbs.tsv           # marked CBSs with CTCF ChIP-seq peaks          
        ├── train.npz                # train set (chr1-16)
        ├── valid.npz                # valid set (chr17-18)
        ├── test.npz                 # test set (chr19-X)                 
        └── total.npz                # feature data for all CBSs
```     


## 3. Run DeepAnchor
DeepAnchor trains a classifier to distinguish insulator-related CBSs from others. After training, the model will be used to predict DeepAnchor score to show the possibility that each CBS belong to insulator-related CBSs. 

To train the classifier model, run command:
```properties
python DeepAnchor.py  work_dir train
```
This will generate a DeepAnchor model which will be saved in DeepAnchor.model in work_dir.

To predict the DeepAnchor score for all CBSs, run command:

```properties
python DeepAnchor.py  work_dir predict
```

This will generate a file *scored_cbs.tsv* that contain all CBSs and their DeepAnchor score. We need to copy this file to ./data/ folder for downstream analyses.

The data columns of *scored_cbs.tsv* are shown below:
|chrom|start|end|strand|score|anchor_score|
|-----|-----|---|------|-----|------------|

*score*: the score for motif scan.
*anchor_score*: the score predicted by DeepAnchor model.

## 4. Run LoopAnchor to make loop prediction

To use LoopAnchor for loop prediction, you should prepare input data arranged as follow:

```
work_dir
    └── raw                   
        └── CTCF_peak.bed.gz
        └── scored_cbs.tsv      
```

Run command:
```properties
python run_LoopAnchor_denovo.py  work_dir
```
In work_dir/LoopAnchor folder, you can find the result LoopAnchor_pred.bedpe which contains all the loops predicted by LoopAnchor.LoopAnchor files is arranged in bedpe format and the last column is the predicted loop intensity.

|chrom1|start1|end1|chrom2|start2|end2|name|score|strand1|strand2|LoopAnchor|
|------|------|----|------|------|----|----|-----|-------|-------|----------|





Here is a complete example. The data can be found in ./data/ folder, but you still need to download some files as shown before.
```properties
python DeepAnchor_input.py ./data/GM12878
python DeepAnchor.py ./data/GM12878 train
python DeepAnchor.py ./data/GM12878 predict
cp ./data/GM12878/scored_cbs.tsv ./data/K562/raw
python run_LoopAnchor_denovo.py ./data/K562
```

## Other code resources

Considering that github may be not available in some regions, we also provide [ReadTheDocs](https://loopanchor.readthedocs.io/en/latest/index.html) and [Bitbuket](https://bitbucket.org/xuhang01/loopanchor/src/main/) resource for LoopAnchor method.


## Landscape availability
We collected 764 available CTCF ChIP-seq data from ENCODE, CistromDB and ChIP-Atlas and use LoopAnchor to predict CTCF-anchored loops. The results are available at UCSC Track Data Hubs (https://genome.ucsc.edu/cgi-bin/hgHubConnect) by entering customized hub URLs https://raw.githubusercontent.com/mulinlab/LoopAnchor/master/loopanchor/data/hubs/hubs_landscape.txt or https://raw.githubusercontent.com/mulinlab/LoopAnchor/master/loopanchor/data/hubs/hubs_all.txt, respectively.


## Citation
Xu, Hang, et al. "Inferring CTCF insulators and anchored loops across human tissues and cell types." bioRxiv (2022).
