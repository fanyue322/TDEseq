cell_aggregate<-function(counts,group,stage,cov,size=20)
{
suppressPackageStartupMessages(library("Seurat"))
meta=data.frame(cell=colnames(counts),group=group,stage=stage)
N=length(unique(group))

pseudolabel<-lapply(1:N,function(x){
g=unique(group)[x]
idx=which(group==g)
if(length(idx)>1)
{
idx=sample(idx)
num=floor(length(idx)/size)+1
pseudolabel=rep(paste0("pseudocell",1:num,"_",g),each=size)[1:length(idx)]
}else{
pseudolabel=rep(paste0("pseudocell",1,"_",g),each=1)
}
return(pseudolabel)
})
pseudolabel=do.call(c,pseudolabel)
meta$pseudolabel=pseudolabel
cov$pseudolabel=pseudolabel
pseudocell=unique(pseudolabel)
pseudo_num=length(pseudocell)
pseudo_counts<-lapply(1:pseudo_num,function(x){
idx=which(pseudolabel==pseudocell[x])
if(length(idx)>1)
{
pseudocount=rowSums(counts[,idx])
}else{
pseudocount=counts[,idx]
}
return(pseudocount)
})
pseudo_counts=do.call(cbind,pseudo_counts)

rownames(pseudo_counts)=rownames(counts)
colnames(pseudo_counts)=pseudocell

idx=duplicated(meta$pseudolabel)
pseudo_meta=meta[which(idx==FALSE),]

for(i in 1:(ncol(cov)-1))
{
for(j in unique(cov$pseudolabel))
{
idx=which(cov$pseudolabel==j)
if(sum(duplicated(cov[idx,i]))!=(length(idx)-1))
{
cov[idx[1],i]=sum(cov[idx,i])
}
}
}

idx=match(pseudo_meta$pseudolabel,colnames(pseudo_counts))		
pseudo_counts=pseudo_counts[,idx]		
seu<-CreateSeuratObject(pseudo_counts)
seu<-NormalizeData(seu)	
seu@meta.data=cbind(seu@meta.data,pseudo_meta)	   
return(seu)
}

logFC_filtering<-function(data,stage)
{
maxFC=rep(0,nrow(data))
stage_idx=sort(unique(stage))
mat=matrix(0,nrow=nrow(data),ncol=length(stage_idx))
for(i in stage_idx)
{
mat[,i+1]=rowMeans(as.matrix(data[,which(stage==i)]))
}
maxFC=apply(mat,1,max)-apply(mat,1,min)
return(maxFC)
}

shape_bstat_lmm<-function(res_dat)
{
shape=rep('NA',nrow(res_dat))
idx=match(c('b.inc', 'b.dec', 'b.cov', 'b.con'),colnames(res_dat))
bstat=res_dat[,idx]
#idx=which((res_dat$shape!='None'))
idx=1:nrow(res_dat)
shape_index=c("Growth","Recession","Trough","Peak")
for(k in idx)
{
if(!is.nan(bstat[k,1]))
{
tmp=bstat[k,]
shape[k]=shape_index[which(tmp==max(tmp)[1])]
}else{
shape[k]=NA
}
}
return(shape)
}

dmvnorm<-function(y,mu,sigma)
{
f<-sum(dnorm(y,mu,sqrt(sigma),log=TRUE))
return(f)
}

reassign_shape<-function(res_dat)
{
shape=rep(NA,nrow(res_dat))
shape_ind=c("Growth","Recession","Trough","Peak")
idx=match(c("aic.inc","aic.dec","aic.cov","aic.con"),colnames(res_dat))
aic=res_dat[,idx]
for(i in 1:nrow(res_dat))
{
aic_tmp=aic[i,]
idx=which(aic_tmp==min(aic_tmp))[1]
shape[i]=shape_ind[idx]
}
return(shape)
}

AIC<-function(dat,time,res)
{
aic=matrix(NA,nrow=nrow(dat),ncol=4)
k=length(unique(time))
n=length(time)
est.inc=t(do.call(cbind,res[seq(2,length(res),9)]))
est.dec=t(do.call(cbind,res[seq(3,length(res),9)]))
est.cov=t(do.call(cbind,res[seq(4,length(res),9)]))
est.con=t(do.call(cbind,res[seq(5,length(res),9)]))
sig.inc=do.call(c,res[seq(6,length(res),9)])
sig.dec=do.call(c,res[seq(7,length(res),9)])
sig.cov=do.call(c,res[seq(8,length(res),9)])
sig.con=do.call(c,res[seq(9,length(res),9)])
for(i in 1:nrow(dat))
{
f=dat[i,]
loglik<-dmvnorm(f,est.inc[i,],sig.inc[i])
aic[i,1]=2*(k+1)-2*loglik
loglik<-dmvnorm(f,est.dec[i,],sig.dec[i])
aic[i,2]=2*(k+1)-2*loglik
loglik<-dmvnorm(f,est.cov[i,],sig.cov[i])
aic[i,3]=2*(k+2)-2*loglik
loglik<-dmvnorm(f,est.con[i,],sig.con[i])
aic[i,4]=2*(k+2)-2*loglik
}
colnames(aic)=c('aic.inc','aic.dec','aic.cov','aic.con')
return(aic)
}

ChangePoint_detection<-function(res_dat,time,res)
{
N=nrow(res_dat)
ChangePoint=rep(NA,N)
shape=res_dat$pattern

est.cov=t(do.call(cbind,res[seq(4,length(res),9)]))
est.con=t(do.call(cbind,res[seq(5,length(res),9)]))

for(i in 1:N)
{
if(shape[i]=='Trough')
{
f=est.cov[i,]
est_t=c()
for(j in unique(time))
{
est_t=c(est_t,mean(f[which(time==j)]))
}
set_dif<-diff(est_t)
pos=which(set_dif>0)[1]
ChangePoint[i]=pos
if(sum(sign(set_dif))== -(length(time)-1))
{
ChangePoint[i]=length(unique(time))
}
}else if(shape[i]=='Peak'){
f=est.con[i,]
est_t=c()
for(j in unique(time))
{
est_t=c(est_t,mean(f[which(time==j)]))
}
set_dif<-diff(est_t)
pos=which(set_dif<0)[1]
ChangePoint[i]=pos
if(sum(sign(set_dif))== (length(time)-1))
{
ChangePoint[i]=length(unique(time))
}
}
}
return(ChangePoint)
}

p_aggregate<-function(pm)
{
pm=as.matrix(pm)
combined_pvalue=rep(NA,nrow(pm))
pm[which(pm<0)]=0
idx=which(is.na(pm[,1]))
if(length(idx)>0)
{
pvalues=pm[-idx,]
res <- setNames(split(pvalues, seq(nrow(pvalues))), rownames(pvalues))
combined_pvalue[-idx] <- unlist(lapply(res, ComputeACAT))
}else{
pvalues=pm
res <- setNames(split(pvalues, seq(nrow(pvalues))), rownames(pvalues))
combined_pvalue <- unlist(lapply(res, ComputeACAT))
}
return(combined_pvalue)
}

GetModelFits<-function(res_dat,res,stage)
{
sig_gene_idx=which(res_dat$SignificantDE=='Yes')
ModelFits=matrix(0,nrow=length(sig_gene_idx),ncol=length(unique(stage)))
rownames(ModelFits)=res_dat$gene[sig_gene_idx]
colnames(ModelFits)=sort(unique(stage))
pattern=c('Growth','Recession','Trough','Peak')
for(i in sig_gene_idx)
{
for(j in sort(unique(stage)))
{
ModelFits[res_dat$gene[i],j+1]=mean(res[[match(res_dat$pattern[i],pattern)]][which(stage==j),i])
}
}
return(ModelFits)
}


ComputeACAT <- function(Pvals, Weights=NULL){
 #### check if there is NA
 if(sum(is.na(Pvals)) > 0){
  stop("Cannot have NAs in the p-values!")
 }## end fi
 
 #### check if Pvals are between 0 and 1
 if((sum(Pvals<0) + sum(Pvals>1)) > 0){
  stop("P-values must be between 0 and 1!")
 }## end fi
 
 #### check if there are pvals that are either exactly 0 or 1.
 is.zero <- (sum(Pvals==0) >= 1)
 is.one <- (sum(Pvals==1) >= 1)
 
 #if(is.zero && is.one){stop("Cannot have both 0 and 1 p-values!")}## end fi
 if(is.zero && is.one){return(NA)}## end fi
 
 if(is.zero){return(1e-300)}## end fi

 if(is.one){
  ##warning("There are p-values that are exactly 1!")
  return(0.9999)
 }## end fi

 #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
 if(is.null(Weights)){
  Weights <- rep(1/length(Pvals), length(Pvals))
 }else if (length(Weights) != length(Pvals)){
  stop("The length of weights should be the same as that of the p-values")
 }else if (sum(Weights<0) > 0){
  stop("All the weights must be positive!")
 }else{
  Weights <- Weights/sum(Weights)
 }## end fi


 #### check if there are very small non-zero p values
 is.small <- (Pvals < 1e-16)
 if(sum(is.small) == 0){
  cct.stat <- sum(Weights*tan((0.5 - Pvals)*pi))
 }else{
  cct.stat <- sum((Weights[is.small]/Pvals[is.small])/pi)
  cct.stat <- cct.stat + sum(Weights[!is.small]*tan((0.5 - Pvals[!is.small])*pi))
 }## end fi

 #### check if the test statistic is very large.
 if(cct.stat > 1e+15){
  pval <- (1/cct.stat)/pi
 }else{
  pval <- 1 - pcauchy(cct.stat)
 }## end fi
 return(pval)
}## end func

shape_bstat_lmm<-function(res_dat)
{
shape=rep('NA',nrow(res_dat))
idx=match(c('b.inc', 'b.dec', 'b.cov', 'b.con'),colnames(res_dat))
bstat=res_dat[,idx]
#idx=which((res_dat$shape!='None'))
idx=1:nrow(res_dat)
shape_index=c("Growth","Recession","Trough","Peak")
for(k in idx)
{
tmp=bstat[k,]
shape[k]=shape_index[which(tmp==max(tmp)[1])]
}
return(shape)
}

#################functions from package CGAM##################

f_ecl = function(evec, ncl, sz) {
	ecl = list()
	st = 1
	ed = 0
	for (icl in 1:ncl) {
		if (length(sz) > 1) {
			szi = sz[icl]
		} else {szi = sz}
		ed = ed + szi
		ecl[[icl]] = evec[st:ed]
		st = ed + 1
	}
	return (ecl)
}



fth2rm = function(th, szs, ycl, N, xcl, p=2, type='b', xtx=NULL, xtx2=NULL, xmat_face=NULL, ones=NULL) {
    ybar = sapply(ycl, mean)
    y = unlist(ycl)
    num = sum(szs^2*ybar^2/(1+szs*th)^2)
    den = sum(y^2) - sum(th*szs^2*ybar^2/(1+szs*th))
    ncl = length(ycl)
    hmat = matrix(0, p, p)
    xtils = list()
    #if (type == 'b') {
    #    n = szs[1]
    #    hmat = xtx - th/(1+n*th)*xtx2
    #} else {
    ones2 = list()
    for(icl in 1:ncl){
        ni = szs[icl]
        xi = xcl[[icl]]
        xm = xi
        #onevec = 1:ni*0+1
        #onemat = onevec%*%t(onevec)
        onemat = ones[[icl]]
        #ones[[icl]] = onemat*th/(1+ni*th)
        ones2[[icl]] = onemat/(1+ni*th)^2
        rinv = diag(ni) - th/(1+ni*th)*onemat
        hmat = hmat + t(xm) %*% rinv %*% xm
        #xtil = t(onevec)%*%xm
        #xtils[[icl]] = xtil
    }
    #oneMat = as.matrix(bdiag(ones))
    #oneMat2 = as.matrix(bdiag(ones2))
    #hmat = xtx - crossprod(xmat_face, oneMat) %*% xmat_face
    
    #ones=ones2=list()
    #    st = 1
    #    ed = 0
    #    for (icl in 1:ncl) {
    #        sz=szs[icl]
    #        ed = ed + sz
    #        onevec = 1:sz*0+1
    #        onemat = onevec%*%t(onevec)
    #        ones[[icl]] = onemat*th/(1+sz*th)
    #        ones2[[icl]] = onemat/(1+sz*th)^2
    #        st = ed + 1
    #    }
    #    oneMat = as.matrix(bdiag(ones))
    #    oneMat2 = as.matrix(bdiag(ones2))
    #    hmat = xtx - crossprod(xmat_face, oneMat) %*% xmat_face
    #}
    #hinv = solve(hmat)
    lmat = chol(hmat)
    hinv = chol2inv(lmat)
    tr = 0
    #if (type == 'b') {
    #    n = szs[1]
        #onevec = 1:n*0+1
        for (icl in 1:ncl) {
            ni = szs[icl]
            xi = xcl[[icl]]
            onevec = 1:ni*0+1
            xtil = t(onevec)%*%xi
            tr = tr + sum(diag(hinv %*% crossprod(xtil)/(1+ni*th)^2))
            #tr = tr + 1/(1+n*th)^2*sum(diag(hinv %*% crossprod(xtil)))
        }
    #    tr = 1/(1+n*th)^2*sum(diag(hinv %*% xtx2))
    #} else {
        # for(icl in 1:ncl) {
        #    ni = szs[icl]
        #    xtil = xtils[[icl]]
        #    tr = tr + sum(diag(hinv %*% crossprod(xtil)/(1+ni*th)^2))
        #}
        #    xtx2ub = crossprod(xmat_face, oneMat2) %*% xmat_face
        #    tr = sum(diag(hinv %*% xtx2ub))
    #}
    rml = 1/2*tr
    obj = (N-p)/2*num/den - 1/2*sum(szs/(1+szs*th)) + rml
    return (obj)
}

fsig = function(thhat, szs, ycl, ncl, N, edf, D, type='b') {
	ybars = sapply(ycl, mean)
    d = min(1.5*edf, D)
	if (type == 'b') {
		sz = N/ncl
		sig2hat = (sum(unlist(ycl)^2) - sz^2*thhat/(1+sz*thhat) * sum(ybars^2))/(N-d-1)
	} else {
		sig2hat = (sum(unlist(ycl)^2) - sum(thhat*szs^2*ybars^2/(1 + szs*thhat)))/(N-d-1)
	}
	return (sig2hat)
}


LogLSFactor<-function(data, scale_factor = scale.factor,num.core=1)
{
    lib.sizes <- colSums(data)
	lib.sizes <- lib.sizes/mean(lib.sizes)
	data <- data/lib.sizes
    res<-pbmcapply::pbmclapply(1:ncol(data), mc.cores = num.core, function(x){
    tryCatch({suppressWarnings(
    res<-normalizedata(d=data[,x],scale_factor=scale_factor)
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
	  data.norm=do.call(cbind,res)
	  rownames(data.norm)=rownames(data)
	  colnames(data.norm)=colnames(data)
	  return(data.norm)
	}
	
normalizedata<-function(d,scale_factor)
{
d<-log2(d+1)
return(d)
}




FastLMM<-function(y,z,eiglist,szs,ncl)		
{

ytil=c()
gtil=c()
ycl = f_ecl(y, ncl, szs) 
zcl = f_ecl(z,ncl,szs)
K_eta_all=c()
for (icl in 1:ncl) {
	ytil=c(ytil,t(eiglist[[icl]]$vectors)%*%ycl[[icl]])
	gtil=c(gtil,t(eiglist[[icl]]$vectors)%*%zcl[[icl]])
	K_eta_all=c(K_eta_all,eiglist[[icl]]$values)
}
			
h2=c()
loglik=c()
 
        interval=c(0,1)
        res<-optimize(likelihood, interval,y_t=ytil,x_t=gtil,K_eta=K_eta_all,maximum=TRUE)
        loglik=c(loglik,res$objective)
        h2=c(h2,res$maximum)

		idx=which(loglik==max(loglik))
        h2=h2[idx]
        sigma<-sigma_delta(h2,gtil,ytil,K_eta_all)
return(list(h2=h2,sigma=sigma))
}

inv_H<-function(theta,K_eta)
{
  z=diag(1/(K_eta*theta+1-theta))
  return(z)
}


exp_sum<-function(theta,y_t,x_t,K_eta)
{
  z=(y_t-x_t*as.numeric(beta_delta(theta,x_t,y_t,K_eta)))^2/(K_eta+(1-theta)/theta)
  return(sum(z))
}

beta_delta<-function(theta,x_t,y_t,K_eta)
{
  ((t(x_t)*(1/(K_eta*theta+1-theta)))%*%y_t)/((t(x_t)*(1/(K_eta*theta+1-theta)))%*%x_t)
}


sigma_delta<-function(theta,x_t,y_t,K_eta)
{
  n=length(y_t)
  1/n*(t(y_t-x_t*as.numeric(beta_delta(theta,x_t,y_t,K_eta)))*(1/(K_eta*theta+1-theta)))%*%(y_t-x_t*as.numeric(beta_delta(theta,x_t,y_t,K_eta)))
}

likelihood<-function(theta,y_t,x_t,K_eta)
{
  n=length(y_t)
  -0.5*(n*log(2*pi)+sum(log(K_eta+(1-theta)/theta))+n+n*log(1/n*exp_sum(theta,y_t,x_t,K_eta)))                          #slow but simple
}

Generate_Kmatrix<-function(szs)
{
eiglist=list()
for(i in 1:length(szs))
{
K=matrix(1,nrow=szs[i],ncol=szs[i])
eig=eigen(K)
eiglist[[i]]=eig
}

return(eiglist)
}


row.multiple<-function(matrix,value)
{
for(i in 1:ncol(matrix))
{
matrix[,i]=matrix[,i]*value[i]
}
return(matrix)
}

