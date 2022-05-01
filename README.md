# LoopAnchor

CTCF is the most important Transcription Factor (TF) for genomic insulation in vertebrate. But it is also a versatile TF that plays several other roles in transcriptional regulation. Here we present DeepAnchor to provide precise description of the genomic/epigenomic patterns surrounding insulation-related CTCF binding sites (CBSs). Generally, DeepAnchor usees cohesin ChIA-PET data, CTCF ChIP-seq data and CBSs predicted by motif scan as input, to train a classifier to distinguish insulation-related CBSs from others. DeepAnchor then calculates a score ranged from [0, 1] for all CBSs, with a larger value meaning more likely to be a positive CTCF insulator.


## Prepare data
To implement DeepAnchor, one should prepare four types of data:
1. Obtain CBSs by motif scanning.
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
        ├── loop.bedpe              
        ├── CTCF_motif.tsv
        ├── dna_feature.npz      
        └── cadd_feature.npz
```
More explainations of data requirement:

CTCF_peak.bed.gz: 
    (1) the ChIP-seq narrowPeak file for CTCF.
    (2) Columns: chrom,chromStart,chromEnd,name,score,strand,signalValue,pValue,qValue,summit.

loop.bedpe:
    (1) the ChIA-PET loop file for Cohesin (RAD21) from the same cell type.
    (2) Columns: chr1,start1,end1,chr2,start2,end2,name,score,strand1,strand2.

CTCF_motif.tsv:
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


It will generate an output folder named DeepAnchor in work_dir:
```
work_dir
    └── DeepAnchor  
        ├── total_anchors.bed        # all anchors extracted from loop.bedpe
        ├── CTCF_peaks.bed           # intersect CTCF_peak.bed.gz and CTCF_motif.tsv
        ├── marked motif.tsv         # mark CBSs with CTCF ChIP-seq peaks          
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

This will generate a file *scored_motif.tsv* that contain all CBSs and their DeepAnchor score. We need to copy this file to ./data/ folder for downstream analyses.

The data columns of *scored_motif.tsv* are shown below:
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
python run_LoopAnchor_denovo.py ./data/K562
```


## Landscape availability
We collected 764 available CTCF ChIP-seq data from ENCODE, CistromDB and ChIP-Atlas and use LoopAnchor to predict CTCF-anchored loops. The results are available at UCSC Track Data Hubs (https://genome.ucsc.edu/cgi-bin/hgHubConnect) by entering customized hub URLs https://raw.githubusercontent.com/mulinlab/LoopAnchor/master/hubs_landscape.txt or https://raw.githubusercontent.com/mulinlab/LoopAnchor/master/hubs_all.txt, respectively.