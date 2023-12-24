#####################################################################
# Package: TDEseq
# Version: 0.0.1
# Modified: 2023-12-9 10:05:02
# Title :  Detecting temporal gene expression changes in the developmental stages of single-cell RNA sequencing studies 
# Authors: Yue Fan and Shiquan Sun
# Contacts: sqsunsph@xjtu.edu.cn;xafanyue@xjtu.edu.cn
#          Xi'an Jiatong University, Department of Biostatistics
######################################################################

#'
#' Fitting the constrained spline model to perform temporal differential expression analysis for time course scRNA-seq data
#' 
#' @param object Gene expression matrix
#' @param fit.method Either 'lmm' or 'lm' to perform (default = lmm).
#' @param pct The percentage of cells where the gene is detected.
#' @param tde.thr A numeric value to indicate significant level of DE genes (default = 0.05).
#' @param lfc Limit testing to genes which show the maximum on average X-fold difference (log-scale) between any two time points. Default is 0.0.
#' @param max.gcells Maximum cell per group. If max.gcells is smaller than the given number of cells in each group, the down-sampling will be active. 
#' @param min.tcells Minimum number of cells in each time points required. 
#' @param tdeseq.method Performing 'cell' or 'pseudocell' strategy. Default is 'cell'.
#' 
#' 
#' @examples
#' data(exampledata)
#' tde <- CreateTDEseqObject(counts = seurat)
#' res<-tdeseq(object = tde)
#' 
#' @export
#' 

tdeseq.default <- function(object,
                  					tde.method = "cell",
                  					sample.id = NULL,
                  					stage.id = NULL,
                  					fit.model = "lmm",
                  					pct = 0.1,
                  					tde.thr = 0.05,
                  					lfc = 0.0,
                  					max.gcells = Inf,
                  					min.tcells = 3,
                  					num.core = 1, 
                  					verbose = FALSE) {
	


	## reordering data ##
	stage_idx=sort(unique(stage.id))
	stage_origin=stage.id
	
	points_order=0
    for(i in stage_idx)
    {
    stage.id[which(stage.id==i)]=points_order
    points_order=points_order+1
    }

    reorder_idx=order(stage.id)
    stage.id=stage.id[reorder_idx]
    stage_origin=stage_origin[reorder_idx]
    object=object[,reorder_idx]
    if(!is.null(sample.id))
    {
    sample.id=sample.id[reorder_idx]
    } 
	
	## filtering time points ##
	num_cell_per_timepoints=table(stage.id)
	idx=names(num_cell_per_timepoints)[which(num_cell_per_timepoints<min.tcells)]
	if(length(idx)>0)
    {
    for(i in idx)
    {
    object=object[,-which(stage.id==i)]
    if(!is.null(sample.id))
    {
    sample.id=sample.id[-which(stage.id==i)]
    }
    stage.id=stage.id[-which(stage.id==i)]
    }
    }
	
	##  downsampling  ##
	if(is.null(max.gcells))
	{
	max.gcells=Inf
	}
	cell_sel=c()
    for(i in unique(sample.id))
    {
    idx=which(sample.id==i)
    if(length(idx)>max.gcells)
    {
    cell_sel=c(cell_sel,sample(idx)[1:max.gcells])
    }else{
    cell_sel=c(cell_sel,idx)
    }
    }
	
	object=object[,cell_sel]
	sample.id=sample.id[cell_sel]
	stage.id=stage.id[cell_sel]
    
	## filtering features ##
	maxFC<-logFC_filtering(object,stage.id)
	idx=which(maxFC>lfc)
	if(length(idx)>0)
	{
	object=object[idx,]
    maxFC=maxFC[idx]
    }
	
	numZeroCounts=rowSums(object==0)
    gene_pct=numZeroCounts/ncol(object)
    idx=which(gene_pct<(1-pct))
    if(length(idx)>0)
	{
    object=object[idx,]
    maxFC=maxFC[idx]
    }
    
	## number of cells and genes
	num_cell <- ncol(object)
	num_gene <- nrow(object)
	genes.use <- rownames(object)
    
	#################
	cat(paste("## ===== TDEseq INPUT INFORMATION ====## \n"))
	cat(paste("## the model fitting: ", fit.model," model\n"))
	cat(paste("## number of total samples: ", ncol(object),"\n"))
	cat(paste("## number of total features: ", nrow(object),"\n"))
	cat(paste("## number of cores: ", num.core,"\n"))
	cat(paste("## ===== END INFORMATION ==== \n"))
	cat("\n")
	
    
	
	## main functions
	if(tde.method == "cell"){
	##*************************************************##
	##   Performing Temporal DE based on Single-Cell   ##
	##*************************************************##
	  if(verbose) cat("# fitting cell-based TDEseq model ... \n")
	  basis=list()
	  basis[[1]]<-basisfunction(x=stage.id,type=1,knots=unique(stage.id),fit.model=fit.model)
      basis[[2]]<-basisfunction(x=stage.id,type=2,knots=unique(stage.id),fit.model=fit.model)  
	  basis[[3]]<-basisfunction(x=stage.id,type=3,knots=unique(stage.id),fit.model=fit.model)
	  basis[[4]]<-basisfunction(x=stage.id,type=4,knots=unique(stage.id),fit.model=fit.model)
		#=====================================
		#res_vc <- parallel::mclapply(seq_len(num_gene), mc.cores = num.core, function(x){
		res.tdeseq <- pbmcapply::pbmclapply(seq_len(num_gene), mc.cores = num.core, function(x){
		  #for each condition get data as y_data
		  tryCatch({suppressWarnings(
		    res <- TDEseq.cell(data = object[x,],
		                       stage = stage.id,
		                       group = sample.id,
		                       z = 0,
		                       fit.model = fit.model,
		                       basis=basis)
		      )
		    }, warning=function(w){ 
		      print(w); return(res);
		    }, error=function(e){
		      print(e); return(NULL);
		    }, finally={
		      #######
		      return(res)
		    }
		  )
		})## end parallel
	}else if(tde.method == "pseudocell"){
	##*************************************************##
	##   Performing Temporal DE based on PseudoCell    ##
	##*************************************************##
	  if(verbose) cat("# fitting pseudocell-based TDEseq model ... \n")
	  res.tdeseq <- pbmcapply::pbmclapply(seq_len(num_gene), mc.cores = num.core, function(x){
	    #for each condition get data as y_data
	    tryCatch({suppressWarnings(
	      res <- TDEseq.pseudocell(data = object[x,],
	                               stage = stage.id,
	                               group = sample.id,
	                               z = 0,
	                               LMM = fit.model,
	                               pct = pct,
	                               threshold = tde.thr,
	                               logFC_threshold = lfc,
	                               max_cells_per_ident = max.gcells,
	                               min_cells_per_timepoints = min.tcells,
								   num.core=num.core)
	    )
	    }, warning=function(w){ 
	      print(w); return(res);
	    }, error=function(e){
	      print(e); return(NULL);
	    }, finally={
	      #######
	      return(res)
	    }
	    )
	  })## end parallel
	}
	
	
	res.tdeseq=do.call(c,res.tdeseq)
	dfTDEseqResults<-do.call(rbind,res.tdeseq[seq(1,length(res.tdeseq),9)])
	p=p_aggregate(dfTDEseqResults)
	dfTDEseqResults$pvalue=p
	rownames(dfTDEseqResults)=genes.use
	dfTDEseqResults$gene=rownames(dfTDEseqResults)
    dfTDEseqResults$padj=p.adjust(p,method='BY')
	dfTDEseqResults$SignificantDE='No'
    dfTDEseqResults$SignificantDE[which(dfTDEseqResults$padj<tde.thr)]='Yes'
	dfTDEseqResults$pattern='None'
	if(fit.model!='lmm')
	{
    aic=AIC(object,stage.id,res.tdeseq)
    dfTDEseqResults=cbind(dfTDEseqResults,aic)
    shape<-reassign_shape(dfTDEseqResults)
	}else{
	dfTDEseqResults$b.inc=do.call(c,res.tdeseq[seq(6,length(res.tdeseq),9)])
	dfTDEseqResults$b.dec=do.call(c,res.tdeseq[seq(7,length(res.tdeseq),9)])
	dfTDEseqResults$b.cov=do.call(c,res.tdeseq[seq(8,length(res.tdeseq),9)])
	dfTDEseqResults$b.con=do.call(c,res.tdeseq[seq(9,length(res.tdeseq),9)])
	shape<-shape_bstat_lmm(dfTDEseqResults)
	}
    dfTDEseqResults$pattern[which(dfTDEseqResults$SignificantDE=='Yes')]=shape[which(dfTDEseqResults$SignificantDE=='Yes')]
    dfTDEseqResults$logFC=maxFC
    ChangePoint<-ChangePoint_detection(dfTDEseqResults,stage.id,res.tdeseq)
	if((!all(is.na(ChangePoint))))
    {
    ChangePoint<-order(unique(stage_origin))[ChangePoint]
    ChangePoint<-stage_idx[ChangePoint]
    }
    dfTDEseqResults$ChangePoint<-ChangePoint
	
	
	
	return(dfTDEseqResults)
}## end function 



#' @param object An TDEseq object 
#'
#' @rdname TDEseq
#' @param fit.method Either 'lmm' or 'lm' to perform (default = lmm).
#' @param pct The percentage of cells where the gene is detected.
#' @param tde.thr A numeric value to indicate significant level of DE genes (default = 0.05).
#' @param lfc Limit testing to genes which show the maximum on average X-fold difference (log-scale) between any two time points. Default is 0.0.
#' @param max.gcells Maximum cell per group. If max.gcells is smaller than the given number of cells in each group, the down-sampling will be active. 
#' @param min.tcells Minimum number of cells in each time points required. 
#' @param tdeseq.method Performing 'cell' or 'pseudocell' strategy. Default is 'cell'.
#' 
#' @method TDEseq TDEseq
#'
#' @examples
#' 
#' \dontrun{
#' data("pbmc_small")
#' object
#' object <- tdeseq(object = object)
#' }
#' 
tdeseq.Assay <- function(object, 
                        	assay = NULL,
                        	slot.use = "data",
                        	features = NULL,
                        	tde.method = "cell",
                        	tde.param = list(sample.var = "group",
                        	                 stage.var = "stage",
                        	                 fit.model = "lmm",
                        	                 pct = 0.1,
                        	                 tde.thr = 0.05,
                        	                 lfc = 0.0,
                        	                 max.gcells = Inf,
                        	                 min.tcells = 3),
                        	num.core = 1,
                        	verbose = TRUE, ...) {
	
  ## data use
	data.use <- GetAssayData(object = object, slot = slot.use)
	meta.data <- GetAssayData(object = object, slot = "meta.data")
	
	if(tde.param$sample.var %in% colnames(meta.data)){
	  sample.id <- meta.data[,tde.param$sample.var,drop=TRUE]
	}else{
	  stop("The variable 'sample.var' is not in the 'meta.data'!")
	}## end fi
	
	if(tde.param$stage.var %in% colnames(meta.data)){
	  stage.id <- meta.data[,tde.param$stage.var,drop=TRUE]
	}else{
	  stop("The variable 'stage.var' is not in the 'meta.data'!")
	}## end fi
	

	if(length(data.use) == 0){
	  counts <- GetAssayData(object = object, slot = "counts")
	  if(length(counts) == 0){
	    stop(paste0("TDEseq::Please provide the slot ", slot.use, " data before running TDEseq function. \n"))
	  }else{
	    data.use <- NormDataTDEseq(data = counts)
	  }## end fi
	  rm(counts);gc();
		
	}## end fi
	##
	features <- features %||% rownames(x = data.use)
	
	
	## run default TDEseq, assume new.data has been ordered by combined p-values
	new.data <- tdeseq(object = data.use[features, , drop=FALSE],
              				tde.method = tde.method,
              				sample.id = sample.id, 
              				stage.id = stage.id,
              				pct = tde.param$pct,
              				tde.thr = tde.param$tde.thr,
              				lfc = tde.param$lfc,
              				max.gcells = tde.param$max.gcells,
              				min.tcells = tde.param$min.tcells,
              				fit.model = tde.param$fit.model,
              				num.core = num.core,
              				verbose = verbose, ...)
					

	## store the scaled data in the slot
	object <- SetAssayData(object = object, slot = 'tde', new.data = new.data)
	## store top number of features in tde slot
}## end func

#' TDEseq: detecting temporal expression changes in time-course RNA sequencing studies.
#' @param object An TDEseq object 
#' 
#' 
#' @param fit.method Either 'lmm' or 'lm' to perform (default = lmm).
#' @param pct The percentage of cells where the gene is detected.
#' @param tde.thr A numeric value to indicate significant level of DE genes (default = 0.05).
#' @param lfc Limit testing to genes which show the maximum on average X-fold difference (log-scale) between any two time points. Default is 0.0.
#' @param max.gcells Maximum cell per group. If max.gcells is smaller than the given number of cells in each group, the down-sampling will be active. 
#' @param min.tcells Minimum number of cells in each time points required. 
#' @param tdeseq.method Performing 'cell' or 'pseudocell' strategy. Default is 'cell'.
#' 
#' @return object An TDEseq object
#' 
#' @author Yue Fan, Shiquan Sun
#' 
#' @examples
#' data(exampledata)
#' res <- tdeseq(object=data)
#' 
#' @export
#' 
tdeseq.TDEseq <- function(object, 
                            assay = 'RNA',
                            slot.use = "data",
                            features = NULL,
                            tde.method = "cell",
                            tde.param = list(sample.var = "group",
                                             stage.var = "stage",
                                             fit.model = "lmm",
                                             pct = 0.1,
                                             tde.thr = 0.05,
                                             lfc = 0.0,
                                             max.gcells = Inf,
                                             min.tcells = 3),
                            num.core = 1,
                            verbose = TRUE, ...) {
	
	## parallel parameter setting
	if(num.core == 1){
		if(slot(object, name="num.core") > 1) {num.core <- slot(object,name = "num.core")}
	}## end fi
	
	## assays
	assay <- assay %||% DefaultAssay(object = object)
	assay.data <- GetAssay(object = object, assay = assay)

	## run main
	assay.data <- tdeseq(object = assay.data,
					slot.use = slot.use,
					features =  features,
					tde.method = tde.method,
					tde.param = tde.param,
					num.core = num.core,
					verbose = verbose, ... )

	## store back the treated data
	object@assays[[assay]] <- assay.data
	return(object)
}## end func


#' Differential expression analysis for SRT
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname tdeseq
#' @export tdeseq
#'
#' @concept data-access
#'
#' @examples
#' tdeseq(object)
#'
tdeseq <- function(object, ...) {
	UseMethod(generic = "tdeseq", object = object)
}## end func



##############################################################################################
#'
#' Normalize raw data
#'
#' Normalize count data per cell and transform to log scale via library size factor
#'
#' @param data Matrix with the raw count data
#' @param scale.factor Scale the data. Default is 1e4
#' @param verbose Print progress
#'
#' @return Returns a matrix with the normalize and log transformed data
#'
#' @import Matrix
#' @importFrom methods as
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
#' mat
#' mat_norm <- LogLibSF(data = mat)
#' mat_norm
#'
NormDataTDEseq <- function(data, scale.factor = 1e4, verbose = TRUE) {
  if (is.data.frame(x = data)) {
    data <- as.matrix(x = data)
  }## end fi
  
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as(object = data, Class = "dgCMatrix")
  }## end fi
  
  ## call Rcpp function to normalize
  if (verbose) {
    cat("Performing log-library size factor normalization\n", file = stderr())
  }## end fi
  
  norm.data <- LogLSFactor(data, scale_factor = scale.factor)
  colnames(x = norm.data) <- colnames(x = data)
  rownames(x = norm.data) <- rownames(x = data)
  return(norm.data)
}## end func



TDEseq.cell<-function(data,
                      stage,
					  group,
					  z=0,
					  fit.model='lmm',
                      basis=basis)
{
if(fit.model!='lmm')
{
fit.incr<-conspline(data,stage,basis=basis[[1]],type=1,test=TRUE)
fit.decr<-conspline(data,stage,basis=basis[[2]],type=2,test=TRUE)
fit.conv<-conspline(data,stage,basis=basis[[3]],type=3,test=TRUE)
fit.conc<-conspline(data,stage,basis=basis[[4]],type=4,test=TRUE)

res<-data.frame(increasing.pvalue=fit.incr$pvalx,decreasing.pvalue=fit.decr$pvalx,convex.pvalue=fit.conv$pvalx,concave.pvalue=fit.conc$pvalx)
results<-list(res=res,est.inc=fit.incr$muhat,est.dec=fit.decr$muhat,est.cov=fit.conv$muhat,est.con=fit.conc$muhat,sig.inc=fit.incr$sighat,sig.dec=fit.decr$sighat,
sig.cov=fit.conv$sighat,sig.con=fit.conc$sighat)
return(results)
}else{
group_numeric=as.character(group)
cout=1
for(i in unique(group))
{
group_numeric[which(group==i)]=cout
cout=cout+1
}
group=as.numeric(group_numeric)


fit.incr<-conespline_lmm(y=data,x=stage,basis=basis[[1]],group,shape=9)
fit.decr<-conespline_lmm(y=data,x=stage,basis=basis[[2]],group,shape=10)
fit.conv<-conespline_lmm(y=data,x=stage,basis=basis[[3]],group,shape=11)
fit.conc<-conespline_lmm(y=data,x=stage,basis=basis[[4]],group,shape=12)



res<-data.frame(increasing.pvalue=fit.incr$pval,decreasing.pvalue=fit.decr$pval,convex.pvalue=fit.conv$pval,concave.pvalue=fit.conc$pval)
results<-list(res=res,est.inc=fit.incr$muhat,est.dec=fit.decr$muhat,est.cov=fit.conv$muhat,est.con=fit.conc$muhat,bstat.inc=fit.incr$bstat,bstat.dec=fit.decr$bstat,
bstat.cov=fit.conv$bstat,bstat.con=fit.conc$bstat)
return(results)

}
}
















#########################################
#             CODE END                  #
#########################################
