# TDEseq
![IMAGE](https://github.com/fanyue322/fanyue322.github.io/blob/master/workflow.png "TDEseq logo")  
TDEseq is implemented as an open source R package for mQTL detecting genes with temporal dynamic expression patterns in time-series scRNA-seq  transcriptomic studies. TDEseq models the relationship between log normalized data and the corresponding time points through constrained additive mixed model and can detect the DE genes as well as its temporal dynamic pattern simultaneously. 


# Installation
TDEseq is implemented as an R package, which can be installed from either GitHub.

```
library(devtools)
install_github("fanyue322/TDEseq")
```

# Usage
The main function is image. You can find the instructions and an example by '?TDEseq'.

## Example
A toy example for testing purposes only:
```
data(ExampleData)
dat <- ExampleData$dat
time <- ExampleData$time
group <- ExampleData$group
res=TDEseq(dat,time,group=group,LMM=TRUE)
```
An example of the outputs TDEseq produces:
```
data(example_results)
```
## Our group

 <https://sqsun.github.io/>.
