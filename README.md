# TDEseq
![TDEseq](https://github.com/fanyue322/fanyue322.github.io/blob/master/workflow_web.png "TDEseq logo")  
TDEseq is implemented as an open source R package for detecting genes with temporal dynamic expression patterns in time-series scRNA-seq transcriptomic studies. TDEseq primarily builds upon the linear additive mixed model (LAMM) framework, with a random effect term to account for correlated samples in time-resolved or time-course scRNA-seq studies. In this model, we typically introduce the quadratic I-splines and cubic C-splines as basis functions, which facilitate the detection of four potential temporal gene expression patterns, i.e., growth, recession, peak, and trough. This vignette will illustrate some uses and visulization of TDEseq.

# Installation
TDEseq is implemented as an R package, which can be installed from either GitHub.

```
library(devtools)
install_github("fanyue322/TDEseq")
```

# Usage
The main function is TDEseq. You can find the instructions and an example by '?TDEseq'.

## Example
A toy example for testing purposes only:
```
data(ExampleData)
dat <- ExampleData$dat
stage <- ExampleData$time
group <- ExampleData$group
metadata<-data.frame(stage=stage,group=group)
rownames(metadata)=colnames(dat)
```
### Linear model version
```
res=TDEseq(X=dat,meta=metadata,LMM=FALSE)
```
### Linear mixed model version
```
res=TDEseq(X=dat,meta=metadata,LMM=TRUE)
```
An example of the outputs TDEseq produces:
```
data(example_results)
```
A tutorial includes main example codes for mouse liver development analysis can be found [here](https://fanyue322.github.io/TDEseq)
## Our group

 <https://sqsun.github.io/>.
