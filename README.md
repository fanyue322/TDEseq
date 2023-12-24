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
The main function is TDEseq. You can find the instructions and an example by '?tdeseq.default'.

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
We show how to create a TDEseqObject object. We can create a TDEseqObject using the count matrix and meta data. Although we provide normalize function to perform log normalization for raw counts data, we recommend the users provided their own normalized scRNA-seq data. The time points information and sample information (for mixed model only) should be contained in the meta.data.
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
Add the parameter setting of TDEseq. 
```
tde_method = "cell"
tde_param = list(sample.var = "batch",
                 stage.var = "stage",
                 fit.model = "lm",
				 tde.thr = 0.05,
		         num.core=10)
tde <- tdeseq(object = tde, tde.param=tde_param)
```
Users need to specify which column in the meta.data corresponds to sample and time points information by setiing `sample.var` and `stage.var`. We set `fit.model="lm"` to perform linear version of TDEseq. Uesr can perform mixed version of TDEseq by setting `fit.model="lmm"`.

### Other options
User can set other parameters to perform some preprocessing steps. This parameter will do four things:
```
tde_param = list(sample.var = "batch",
                 stage.var = "stage",
                 fit.model = "lm",
                 pct = 0.1,
                 tde.thr = 0.05,
                 lfc = 0.1,
                 max.gcells = Inf,
                 min.tcells = 3,
		 num.core=10)
tde <- tdeseq(object = tde, tde.param=tde_param)
```
1. Remove time points with too few cells by setting `min.tcells`. Here, time points with less than 3 cells will be removed.
2. Filter genes that are only expressed in a few cells by setting `pct`. Here, genes with more than 90% of zero counts will be filtered out.
3. Filter genes that show small average X-fold difference (log-scale) between any two time points by setting `lfc`. Here, we limit testing to genes which show at least 0.1-fold difference between any two time points.
4. Downsample cells by setting `max.gcells`. If max.gcells is smaller than the given number of cells in a sample, the down-sampling will be active. Here, we do not perform downsampling by setting `max.gcells=Inf`.

### Get results
The results of TDEseq analysis are stored in TDEseqObject. User can obtain the results by
```
## Get the results of TDEseq analysis for each gene
result<-GetTDEseqAssayData(tde,'tde')  
```
A tutorial includes main example codes for mouse liver development analysis can be found [here](https://fanyue322.github.io/TDEseq)
## Our group

 <https://sqsun.github.io/>.
