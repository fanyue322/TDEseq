# TDEseq
![TDEseq](https://github.com/fanyue322/fanyue322.github.io/blob/master/workflow_web.png "TDEseq logo")  
TDEseq is implemented as an open source R package for detecting genes with temporal dynamic expression patterns in time-series scRNA-seq transcriptomic studies. TDEseq primarily builds upon the linear additive mixed model (LAMM) framework, with a random effect term to account for correlated samples in time-resolved or time-course scRNA-seq studies. In this model, we typically introduce the quadratic I-splines and cubic C-splines as basis functions, which facilitate the detection of four potential temporal gene expression patterns, i.e., growth, recession, peak, and trough. This vignette will illustrate some uses and visulization of TDEseq.

# Installation
TDEseq is implemented as an R package, which can be installed from GitHub.

```
library(devtools)
install_github("fanyue322/TDEseq")
```

# Usage
The main function is TDEseq. You can find the instructions and an example by '?tdeseq'.

## Example
We demonstrate the use of TDEseq to an example simulated time course scRNA-seq data that are here, which are included in the TDEseq package. This toy example is used for testing purposes only:

### Load the simulated data
```
data(exampledata)
seurat
#> An object of class Seurat 
#> 200 features across 1246 samples within 1 assay 
#> Active assay: RNA (200 features, 0 variable features)
```
### Create a TDEseq object
We show how to create a TDEseqObject object. We can create a TDEseqObject using the count matrix and meta data. Although we provide normalize function to perform log normalization for raw counts data, we recommend the users provided their own normalized scRNA-seq data. 
```
counts=Seurat::GetAssayData(seurat,'counts')
norm.data=Seurat::GetAssayData(seurat,'data')
tde <- CreateTDEseqObject(counts = counts, data=norm.data, meta.data = seurat@meta.data)
```
Note: the time points and sample information must be contained in the meta data.

Alternatively, TDEseqObject can be created directly from a Seurat object 
```
tde <- CreateTDEseqObject(counts = seurat)
```
or a sce object
```
sce <- SingleCellExperiment::SingleCellExperiment(list(counts=counts,logcounts=data.norm),
                            colData=data.frame(label=colnames(counts)),
                            rowData=data.frame(length=rownames(counts)),
                            metadata=list(study="GSE111111"))
tde <- CreateTDEseqObject(counts = sce)
```
### Fit TDEseq using simulated data
```
res=TDEseq(X=dat,meta=metadata,LMM=FALSE)
```
### Linear mixed model version
```
res=TDEseq(X=dat,meta=metadata,LMM=TRUE)
```
A tutorial includes main example codes for mouse liver development analysis can be found [here](https://fanyue322.github.io/TDEseq)
## Our group

 <https://sqsun.github.io/>.
