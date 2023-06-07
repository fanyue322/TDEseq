########################################################################################################################
# Package: TDEseq
# Version: 1.0
# Date   : 2022-11-22
# Title  : TDEseq: detecting temporal gene expression changes in the developmental stages of single-cell RNA sequencing studies
# Authors: Yue Fan and Shiquan Sun
# Contact: xafanyue@xjtu.edu.cn
#          Xian Jiaotong university, Department of Public Health
########################################################################################################################

#' TDEseq: detecting temporal gene expression changes in the developmental stages of single-cell RNA sequencing studies.
#' 
#' @param data The log normalization of expression values, with genes in rows and cells in columns.
#' @param stage A vector of stage index or time points.
#' @param group A vector to represent the random effect.
#' @param LMM A bool value to indicate whether to perform LMM or LM analysis (default = FALSE).
#' @param threshold A numeric value to indicate significant level of DE genes (default = 0.05).
#' @return \item{gene}{The name of genes}
#' @return \item{pval}{P value for the fixed effect}
#' @return \item{padj}{Adjusted P value for the fixed effect}
#' @return \item{pattern}{The most appropriate pattern for each DE genes}
#' @return \item{ChangePoint}{The change point for DE genes, only exists for peak or recession DE genes}
#' @author Yue Fan, Shiquan Sun
#' @examples
#' data(exampledata)
#' data=seurat@assays$RNA@data
#' stage=seurat@meta.data$stage
#' res<-TDEseq(data,stage)
#' @export

TDEseq<-function(data,stage,z=0,group=NULL,verbose=TRUE,LMM=FALSE,threshold=0.05) {

numCell=ncol(data)
numVar=nrow(data)
cat(paste("## number of total genes: ", numVar,"\n"))
cat(paste("## number of total cells: ", numCell,"\n"))

if(numCell!=length(stage))
{
cat(paste("## Error: the stage information and the cells are not matched!"))
}

if(is.null(rownames(data)))
{
cat(paste("## Error: Please provide the gene names as the row names of the gene expression matrix!"))
}

if(LMM==FALSE)
{

res<-TDEseq_lm(data,stage,z=z,verbose)
res_dat=res$res
p=p_aggregate(res_dat[,2:5])
res_dat$pval=p
res_dat$padj=p.adjust(p,method='fdr')
aic=AIC(data,stage,res)
res_dat=cbind(res_dat,aic)
res_dat$SignificantDE='No'
res_dat$SignificantDE[which(res_dat$padj<threshold)]='Yes'
shape<-reassign_shape(res_dat)
res_dat$pattern=shape
ChangePoint<-ChangePoint_detection(res_dat,stage,res)
ChangePoint<-order(unique(stage))[ChangePoint]
res_dat$ChangePoint<-ChangePoint


}else if(LMM==TRUE){

res<-TDEseq_lmm(data,stage,group,z=z,verbose)
res_dat<-res$res
p<-p_aggregate(res_dat[,2:5])
res_dat$pval<-p
res_dat$padj<-p.adjust(p,method='fdr')
res_dat$SignificantDE='No'
res_dat$SignificantDE[which(res_dat$padj<threshold)]='Yes'
bstat<-res$bstat
colnames(bstat)<-c('b.inc', 'b.dec', 'b.cov', 'b.con')
res_dat<-cbind(res_dat,bstat)
shape<-shape_bstat_lmm(res_dat)
res_dat$pattern=shape
ChangePoint<-ChangePoint_detection(res_dat,stage,res)
ChangePoint<-order(unique(stage))[ChangePoint]
res_dat$ChangePoint=ChangePoint
}

return(res_dat)

}


TDEseq_lm<-function(dat,stage,z=0,verbose=TRUE)
{
est.inc=matrix(0,nrow=0,ncol=length(stage))
est.dec=matrix(0,nrow=0,ncol=length(stage))
est.con=matrix(0,nrow=0,ncol=length(stage))
est.cov=matrix(0,nrow=0,ncol=length(stage))
sig_est.inc=c()
sig_est.dec=c()
sig_est.con=c()
sig_est.cov=c()
res_dat=data.frame()
x=stage
for(iVar in 1:nrow(dat))
{
  if(verbose==TRUE)
  {
  print(paste0('processing: ',rownames(dat)[iVar]))
  }
  N=ncol(dat)
  Expr=dat[iVar,] 
  
  y<-Expr
 
  fit.incr<-conspline(y,x,type=1,knots=unique(x),test=TRUE,nsim=100,zmat=z)
  est.inc=rbind(est.inc,as.numeric(fit.incr$muhat))
  sig_est.inc=c(sig_est.inc,fit.incr$sighat)
 fit.decr<-conspline(y,x,type=2,knots=unique(x),test=TRUE,nsim=100,zmat=z)
  est.dec=rbind(est.dec,as.numeric(fit.decr$muhat))
  sig_est.dec=c(sig_est.dec,fit.decr$sighat)
 fit.conv<-conspline(y,x,type=3,knots=unique(x),test=TRUE,nsim=100,zmat=z)
  est.cov=rbind(est.cov,as.numeric(fit.conv$muhat))
  sig_est.cov=c(sig_est.cov,fit.conv$sighat)
 fit.conc<-conspline(y,x,type=4,knots=unique(x),test=TRUE,nsim=100,zmat=z)
  est.con=rbind(est.con,as.numeric(fit.conc$muhat))
  sig_est.con=c(sig_est.con,fit.conc$sighat)

res<-data.frame(gene=rownames(dat)[iVar],increasing.pvalue=fit.incr$pvalx,
decreasing.pvalue=fit.decr$pvalx,convex.pvalue=fit.conv$pvalx,concave.pvalue=fit.conc$pvalx)
  res_dat=rbind(res_dat,res)
}

results<-list(res=res_dat,est.inc=est.inc,est.dec=est.dec,est.cov=est.cov,est.con=est.con,sig.inc=sig_est.inc,sig.dec=sig_est.dec,
sig.con=sig_est.con,sig.cov=sig_est.cov)
return(results)
}


TDEseq_lmm<-function(data,stage,group,z=0,verbose=TRUE)
{
est.inc=matrix(0,nrow=0,ncol=length(stage))
est.dec=matrix(0,nrow=0,ncol=length(stage))
est.con=matrix(0,nrow=0,ncol=length(stage))
est.cov=matrix(0,nrow=0,ncol=length(stage))
sig_est.inc=c()
sig_est.dec=c()
sig_est.con=c()
sig_est.cov=c()
bstat=matrix(0,nrow=nrow(dat),ncol=4)
res_dat=data.frame()
x=stage
for(iVar in 1:nrow(dat))
{
  if(verbose==TRUE)
  {
  print(paste0('processing: ',rownames(dat)[iVar]))
  }
   N=ncol(dat)
  Expr=dat[iVar,] 
  
  y<-Expr

if(length(z==1))
{
fit.incr<-conespline_lmm(x,y,group,shape=9)
est.inc=rbind(est.inc,as.numeric(fit.incr$muhat))
sig_est.inc=c(sig_est.inc,fit.incr$sig2hat)
bstat[iVar,1]=fit.incr$bstat
fit.decr<-conespline_lmm(x,y,group,shape=10)
est.dec=rbind(est.dec,as.numeric(fit.decr$muhat))
sig_est.dec=c(sig_est.dec,fit.decr$sig2hat)
bstat[iVar,2]=fit.decr$bstat
fit.conv<-conespline_lmm(x,y,group,shape=11)
est.cov=rbind(est.cov,as.numeric(fit.conv$muhat))
sig_est.cov=c(sig_est.cov,fit.conv$sig2hat)
bstat[iVar,3]=fit.conv$bstat
fit.conc<-conespline_lmm(x,y,group,shape=12)
est.con=rbind(est.con,as.numeric(fit.conc$muhat))
sig_est.con=c(sig_est.con,fit.conc$sig2hat)
bstat[iVar,4]=fit.conc$bstat
}else{
fit.incr<-conespline_lmm_cov(x,y,z=z,group,shape=9)
est.inc=rbind(est.inc,as.numeric(fit.incr$muhat))
sig_est.inc=c(sig_est.inc,fit.incr$sig2hat)
bstat[iVar,1]=fit.incr$bstat
fit.decr<-conespline_lmm_cov(x,y,z=z,group,shape=10)
est.dec=rbind(est.dec,as.numeric(fit.decr$muhat))
sig_est.dec=c(sig_est.dec,fit.decr$sig2hat)
bstat[iVar,2]=fit.decr$bstat
fit.conv<-conespline_lmm_cov(x,y,z=z,group,shape=11)
est.cov=rbind(est.cov,as.numeric(fit.conv$muhat))
sig_est.cov=c(sig_est.cov,fit.conv$sig2hat)
bstat[iVar,3]=fit.conv$bstat
fit.conc<-conespline_lmm_cov(x,y,z=z,group,shape=12)
est.con=rbind(est.con,as.numeric(fit.conc$muhat))
sig_est.con=c(sig_est.con,fit.conc$sig2hat)
bstat[iVar,4]=fit.conc$bstat
}
res<-data.frame(gene=rownames(dat)[iVar],increasing.pvalue=fit.incr$pval,
decreasing.pvalue=fit.decr$pval,convex.pvalue=fit.conv$pval,concave.pvalue=fit.conc$pval)
res_dat=rbind(res_dat,res)

}
results<-list(res=res_dat,est.inc=est.inc,est.dec=est.dec,est.cov=est.cov,est.con=est.con,sig.inc=sig_est.inc,sig.dec=sig_est.dec,
sig.con=sig_est.con,sig.cov=sig_est.cov,bstat=bstat)
return(results)
}

p_aggregate<-function(pm)
{
pm=as.matrix(pm)
pm[which(pm<0)]=0
idx=which(is.na(pm[,1]))
if(length(idx)>0)
{
pvalues=pm[-idx,]
}else{
pvalues=pm
}
res <- setNames(split(pvalues, seq(nrow(pvalues))), rownames(pvalues))
combined_pvalue <- unlist(lapply(res, ComputeACAT))
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

###################LM function###############################

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
for(i in 1:nrow(dat))
{
f=dat[i,]
loglik<-dmvnorm(f,res$est.inc[i,],res$sig.inc[i])
aic[i,1]=2*(k+1)-2*loglik
loglik<-dmvnorm(f,res$est.dec[i,],res$sig.dec[i])
aic[i,2]=2*(k+1)-2*loglik
loglik<-dmvnorm(f,res$est.cov[i,],res$sig.cov[i])
aic[i,3]=2*(k+2)-2*loglik
loglik<-dmvnorm(f,res$est.con[i,],res$sig.con[i])
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
for(i in 1:N)
{
if(shape[i]=='Trough')
{
f=res$est.cov[i,]
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
f=res$est.con[i,]
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

############LMM function########################
conespline_lmm<-function(x,y,group,shape=9,test=TRUE,nsim=100)
{ ## adjust
xmat=as.matrix(x)
n = length(y)
id=group
    szs = unname(table(id)) ### numbers of individuals in each group
    #print (id)
    ncl = length(szs) ###numbers of groups
	balanced = FALSE
	ycl = f_ecl(y, ncl, szs)  #observations for each group
	sm = 1e-7 
	capl = length(xmat) / n  ##numbers of predictors
	delta = NULL
	varlist = NULL
	xid1 = NULL; xid2 = NULL; xpos2 = 0  ##position of nonlinear term
	knotsuse = list(); numknotsuse = NULL
	mslst = list()
#new:
    knots=list()
	capm = 0
	capms = 0
	numknots=0
	knots[[1]]=0
	space="E"
	shapes=shape
	    del1_ans = makedelta(xmat[, 1], shape, numknots[1], knots[[1]], space = space[1])
		del1 = del1_ans$amat
		knotsuse[[1]] = del1_ans$knots 
		mslst[[1]] = del1_ans$ms
		numknotsuse = c(numknotsuse, length(del1_ans$knots))
        m1 = length(del1) / n
        var1 = 1:m1*0 + 1
		xpos1 = xpos2 + 1
		xpos2 = xpos2 + m1
		xid1 = c(xid1, xpos1)
		xid2 = c(xid2, xpos2)
		delta = del1
        varlist = var1
		
		 xvec = NULL
		 
if(shape==9 | shape==10)
{
bigmat = rbind(1:n*0 + 1, delta)
np = 1 + capms
}else{
bigmat <- rbind(1:n*0 + 1, t(xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13]), delta)
np <- 1 + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) + capms
}
capm <- length(delta) / n - capms
		
		#new: capm is the number of columns of edges for constrained x's
		capm = length(delta) / n - capms
		
		zvec = y
        gmat = t(bigmat)
		
	dsend = gmat[, (np + 1):(np + capm), drop = FALSE]
        zsend = gmat[, 1:np, drop = FALSE]
		 ans = coneB(zvec, dsend, zsend)
            edf = ans$df
            face = ans$face
            bh = coef(ans)
			    if (any(round(bh[1:np],6) < 0)) {
                pos = (1:np)[which(round(bh[1:np],6) < 0)]
                face = unique(c(pos, face))
            }
			
        dd = t(bigmat[face, ,drop = FALSE])
		
		xms = ones = list()
        st = 1
        ed = 0
		
		for (icl in 1:ncl) {
            sz = szs[icl]
            ed = ed + sz
            xms[[icl]] = dd[st:ed, ,drop=F]
            onevec = 1:sz*0+1
            onemat = onevec%*%t(onevec)
            ones[[icl]] = onemat
            st = ed + 1
        }
		
		muhat = t(bigmat) %*% bh
		oldmu = muhat
#########update mu and sigma iterately##########		
		diff = 10
		nrep = 0
		while (diff > 1e-7 & nrep < 10) {
		
		nrep = nrep + 1
		evec = y - muhat    ##residuals a+e
		ecl = f_ecl(evec, ncl, szs)  ## residuals by group
        mod.lmer = NULL
		###estimate variance of a by lmer package### 
	#	mod.lmer = lmer(evec~-1+(1|id), REML=reml)
    #    thhat = summary(mod.lmer)$optinfo$val^2
		###estimate variance by grid search###
		ansi = try(ansi0<-uniroot(fth2rm, c(1e-10, 1e+3), szs=szs, ycl=ecl, N=n, xcl=xms, p=edf, type='ub', xtx=xtx, xtx2=xtx2, xmat_face=dd, ones=ones), silent=TRUE)
        if (class(ansi) == "try-error") {
            thhat = 0
        } else {
            thhat = ansi$root
        }
		
		type = "ub"
############update mu gaven a ############		
	ytil = NULL 
#gtil is edges
			gtil = NULL
			st = 1
			ed = 0
			sz = max(szs)
            pos = which(szs == sz)[1]
            oneMat = ones[[pos]]
			vi = diag(sz) + oneMat*thhat  ##covariance matrix
            covi = vi
            umat = t(chol(covi))
            uinv = solve(umat)
            #uinv0 is used for unbalanced
            uinv0 = uinv
			######L^-1*y=L^-1*(mu+xb)+e  e~N(0,I) #########
			for (icl in 1:ncl) {
				sz = szs[icl]
                uinv = uinv0[1:sz, 1:sz, drop=FALSE]
				yi = ycl[[icl]]
				ytil = c(ytil, uinv %*% as.matrix(yi, ncol=1))
				ed = ed + sz
				gtil = rbind(gtil, uinv %*% gmat[st:ed, ,drop=F])
				st = ed + 1
			}
			#####weighted coneB #########
			dsend = gtil[, (np + 1):(np + capm), drop = FALSE]
            zsend = gtil[, 1:np, drop = FALSE]
            ans = coneB(ytil, dsend, vmat = zsend, face=face)
            edf = ans$df
            face = ans$face
            bh = coef(ans)
			    if (any(round(bh[1:np],6) < 0)) {
                pos = (1:np)[which(round(bh[1:np],6) < 0)]
                face = unique(c(pos, face))
                }
			muhat = t(bigmat) %*% bh
			diff = mean((oldmu - muhat)^2)
			oldmu = muhat
            dd = t(bigmat[face, ,drop = FALSE])
            dd2 = gtil[,face,drop=FALSE]
                xms = list()
                st = 1
                ed = 0
                for (icl in 1:ncl) {
                    sz = szs[icl]
                    ed = ed + sz
                    xms[[icl]] = dd[st:ed, ,drop=F]
                    st = ed + 1
                }
			}
		ebars = sapply(ecl, mean)
		sig2hat = fsig(thhat, szs, ecl, ncl, N=n, edf=edf, D=nrow(bigmat), type=type)
		siga2hat = sig2hat * thhat 
		ahat = ebars*szs*thhat/(1+szs*thhat)
	#################testing####################
	if(test)
	{
	        ytil=NULL
	        gtil = NULL
			st = 1
			ed = 0
			sz = max(szs)
            pos = which(szs == sz)[1]
            oneMat = ones[[pos]]
			vi = diag(sz) + oneMat*thhat  ##covariance matrix
            covi = vi
            umat = t(chol(covi))
            uinv = solve(umat)
            #uinv0 is used for unbalanced
            uinv0 = uinv
			######L^-1*y=L^-1*(mu+xb)+e  e~N(0,I) #########
			for (icl in 1:ncl) {
				sz = szs[icl]
                uinv = uinv0[1:sz, 1:sz, drop=FALSE]
				yi = ycl[[icl]]
				ytil = c(ytil, uinv %*% as.matrix(yi, ncol=1))
				ed = ed + sz
				gtil = rbind(gtil, uinv %*% gmat[st:ed, ,drop=F])
				st = ed + 1
			}
			#####weighted coneB #########
			dsend = gtil[, (np + 1):(np + capm), drop = FALSE]
            zsend = gtil[, 1:np, drop = FALSE]
			
			yhat=gtil%*%bh                                                    #34.20664
			pmat=zsend%*%solve(t(zsend)%*%zsend)%*%t(zsend)
			th0=pmat%*%ytil
			sse0=sum((ytil-th0)^2)
			sse1=sum((ytil-yhat)^2)
		bstat=(sse0-sse1)/sse0
		m=ncol(gtil)
		mdist=1:(m+1)*0
		k0=dim(zsend)[2]
		for(isim in 1:nsim){
			ysim=rnorm(n)
			asim=coneB(ysim,dsend,zsend)
			df0=asim$df-k0
			mdist[df0+1]=mdist[df0+1]+1
		}
		mdist=mdist/nsim
		ps=mdist[1]
		for(d in 1:m){
			ps=ps+pbeta(bstat,d/2,(n-d-k0)/2)*mdist[d+1]
		}
		pval=1-ps
    }
	rslt = list(muhat = muhat,  bh = bh,ahat = ahat, sig2hat = sig2hat, siga2hat = siga2hat, thhat = thhat, bigmat = bigmat,pval=pval,bstat=bstat)
    return (rslt)
}

conespline_lmm_cov<-function(x,y,z,group,shape=9,test=TRUE,nsim=100)
{ ## adjust
xmat=as.matrix(x)
n = length(y)
zmat=as.matrix(z)
id=group
    szs = unname(table(id)) ### numbers of individuals in each group
    #print (id)
    ncl = length(szs) ###numbers of groups
	balanced = FALSE
	ycl = f_ecl(y, ncl, szs)  #observations for each group
	sm = 1e-7 
	capl = length(xmat) / n  ##numbers of predictors
	capk = length(zmat) / n  ##numbers of covariates
	delta = NULL
	varlist = NULL
	xid1 = NULL; xid2 = NULL; xpos2 = 0  ##position of nonlinear term
	knotsuse = list(); numknotsuse = NULL
	mslst = list()
#new:
    knots=list()
	capm = 0
	capms = 0
	numknots=0
	knots[[1]]=0
	space="E"
	shapes=shape
	    del1_ans = makedelta(xmat[, 1], shape, numknots[1], knots[[1]], space = space[1])
		del1 = del1_ans$amat
		knotsuse[[1]] = del1_ans$knots 
		mslst[[1]] = del1_ans$ms
		numknotsuse = c(numknotsuse, length(del1_ans$knots))
        m1 = length(del1) / n   ####number of nonlinear terms
        var1 = 1:m1*0 + 1
		xpos1 = xpos2 + 1
		xpos2 = xpos2 + m1
		xid1 = c(xid1, xpos1)
		xid2 = c(xid2, xpos2)
		delta = del1
        varlist = var1
		
		 xvec = NULL
		 
if(shape==9 | shape==10)
{
# bigmat = rbind(1:n*0 + 1, delta)
# np = 1 + capms
bigmat = rbind(1:n*0 + 1, t(zmat), delta)
np = 1 + capk + capms
}else{
# bigmat <- rbind(1:n*0 + 1, t(xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13]), delta)
# np <- 1 + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) + capms
xvec = t(xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13])
bigmat = rbind(1:n*0 + 1, t(zmat), xvec, delta)
np = 1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13)  + capms
}

		
		#new: capm is the number of columns of edges for constrained x's
		capm = length(delta) / n - capms
		
		zvec = y
        gmat = t(bigmat)
		
	dsend = gmat[, (np + 1):(np + capm), drop = FALSE]
        zsend = gmat[, 1:np, drop = FALSE]
		 ans = coneB(zvec, dsend, zsend)
            edf = ans$df
            face = ans$face
            bh = coef(ans)
			    if (any(round(bh[1:np],6) < 0)) {
                pos = (1:np)[which(round(bh[1:np],6) < 0)]
                face = unique(c(pos, face))
            }
			
        dd = t(bigmat[face, ,drop = FALSE])
		
		xms = ones = list()
        st = 1
        ed = 0
		
		for (icl in 1:ncl) {
            sz = szs[icl]
            ed = ed + sz
            xms[[icl]] = dd[st:ed, ,drop=F]
            onevec = 1:sz*0+1
            onemat = onevec%*%t(onevec)
            ones[[icl]] = onemat
            st = ed + 1
        }
		
		muhat = t(bigmat) %*% bh
		oldmu = muhat
#########update mu and sigma iterately##########		
		diff = 10
		nrep = 0
		while (diff > 1e-7 & nrep < 10) {
		
		nrep = nrep + 1
		evec = y - muhat    ##residuals a+e
		ecl = f_ecl(evec, ncl, szs)  ## residuals by group
        mod.lmer = NULL
		###estimate variance of a by lmer package### 
	#	mod.lmer = lmer(evec~-1+(1|id), REML=reml)
    #    thhat = summary(mod.lmer)$optinfo$val^2
		###estimate variance by grid search###
		ansi = try(ansi0<-uniroot(fth2rm, c(1e-10, 1e+3), szs=szs, ycl=ecl, N=n, xcl=xms, p=edf, type='ub', xtx=xtx, xtx2=xtx2, xmat_face=dd, ones=ones), silent=TRUE)
        if (class(ansi) == "try-error") {
            thhat = 0
        } else {
            thhat = ansi$root
        }
		
		type = "ub"
############update mu gaven a ############		
	        ytil = NULL 
            #gtil is edges
			gtil = NULL
			st = 1
			ed = 0
			sz = max(szs)
            pos = which(szs == sz)[1]
            oneMat = ones[[pos]]
			vi = diag(sz) + oneMat*thhat  ##covariance matrix
			# eig=eigen(vi)
			# emat=(eig$vectors%*%sqrt(diag(eig$values)))
            # uinv=solve(emat)
            covi = vi
            umat = t(chol(covi))
            uinv = solve(umat)
            #uinv0 is used for unbalanced
            uinv0 = uinv
			######L^-1*y=L^-1*(mu+xb)+e  e~N(0,I) #########
			for (icl in 1:ncl) {
				sz = szs[icl]
                uinv = uinv0[1:sz, 1:sz, drop=FALSE]
				yi = ycl[[icl]]
				ytil = c(ytil, uinv %*% as.matrix(yi, ncol=1))
				ed = ed + sz
				gtil = rbind(gtil, uinv %*% gmat[st:ed, ,drop=F])
				st = ed + 1
			}
			#####weighted coneB #########
			dsend = gtil[, (np + 1):(np + capm), drop = FALSE]
            zsend = gtil[, 1:np, drop = FALSE]
            ans = coneB(ytil, dsend, vmat = zsend, face=face)
            edf = ans$df
            face = ans$face
            bh = coef(ans)
			    if (any(round(bh[1:np],6) < 0)) {
                pos = (1:np)[which(round(bh[1:np],6) < 0)]
                face = unique(c(pos, face))
                }
			muhat = t(bigmat) %*% bh
			diff = mean((oldmu - muhat)^2)
			oldmu = muhat
            dd = t(bigmat[face, ,drop = FALSE])
            dd2 = gtil[,face,drop=FALSE]
                xms = list()
                st = 1
                ed = 0
                for (icl in 1:ncl) {
                    sz = szs[icl]
                    ed = ed + sz
                    xms[[icl]] = dd[st:ed, ,drop=F]
                    st = ed + 1
                }
			}
		ebars = sapply(ecl, mean)
		sig2hat = fsig(thhat, szs, ecl, ncl, N=n, edf=edf, D=nrow(bigmat), type=type)
		siga2hat = sig2hat * thhat 
		ahat = ebars*szs*thhat/(1+szs*thhat)
	#################testing####################
	if(test)
	{
	        ytil=NULL
	        gtil = NULL
			st = 1
			ed = 0
			sz = max(szs)
            pos = which(szs == sz)[1]
            oneMat = ones[[pos]]
			vi = diag(sz) + oneMat*thhat  ##covariance matrix
			covi = vi
            umat = t(chol(covi))
            uinv = solve(umat)
            # eig=eigen(vi)
			# emat=(eig$vectors%*%sqrt(diag(eig$values)))
            # uinv=solve(emat)
            #uinv0 is used for unbalanced
            uinv0 = uinv
			######L^-1*y=L^-1*(mu+xb)+e  e~N(0,I) #########
			for (icl in 1:ncl) {
				sz = szs[icl]
                uinv = uinv0[1:sz, 1:sz, drop=FALSE]
				yi = ycl[[icl]]
				ytil = c(ytil, uinv %*% as.matrix(yi, ncol=1))
				ed = ed + sz
				gtil = rbind(gtil, uinv %*% gmat[st:ed, ,drop=F])
				st = ed + 1
			}
			#####weighted coneB #########
			dsend = gtil[, (np + 1):(np + capm), drop = FALSE]
            zsend = gtil[, 1:np, drop = FALSE]
			
			yhat=gtil%*%bh                                                    
			pmat=zsend%*%solve(t(zsend)%*%zsend)%*%t(zsend)
			th0=pmat%*%ytil
			sse0=sum((ytil-th0)^2)
			sse1=sum((ytil-yhat)^2)
		bstat=(sse0-sse1)/sse0
		m=ncol(gtil)
		mdist=1:(m+1)*0
		k0=dim(zsend)[2]
		for(isim in 1:nsim){
			ysim=rnorm(n)
			asim=coneB(ysim,dsend,zsend)
			df0=asim$df-k0
			mdist[df0+1]=mdist[df0+1]+1
		}
		mdist=mdist/nsim
		ps=mdist[1]
		for(d in 1:m){
			ps=ps+pbeta(bstat,d/2,(n-d-k0)/2)*mdist[d+1]
		}
		pval=1-ps
    }
	rslt = list(muhat = muhat,  bh = bh,ahat = ahat, sig2hat = sig2hat, siga2hat = siga2hat, thhat = thhat, bigmat = bigmat,pval=pval,bstat=bstat)
    return (rslt)
}

coneB_lmm<-function(gmat,zvec,bigmat,np,capm,szs,ncl,ycl)
{
 dsend = gmat[, (np + 1):(np + capm), drop = FALSE]
        zsend = gmat[, 1:np, drop = FALSE]
		 ans = coneB(zvec, dsend, zsend)
            edf = ans$df
            face = ans$face
            bh = coef(ans)
			    if (any(round(bh[1:np],6) < 0)) {
                pos = (1:np)[which(round(bh[1:np],6) < 0)]
                face = unique(c(pos, face))
            }
			
        dd = t(bigmat[face, ,drop = FALSE])
		
		xms = ones = list()
        st = 1
        ed = 0
		
		for (icl in 1:ncl) {
            sz = szs[icl]
            ed = ed + sz
            xms[[icl]] = dd[st:ed, ,drop=F]
            onevec = 1:sz*0+1
            onemat = onevec%*%t(onevec)
            ones[[icl]] = onemat
            st = ed + 1
        }
		
		muhat = t(bigmat) %*% bh
		oldmu = muhat
#########update mu and sigma iterately##########		
		diff = 10
		nrep = 0
		while (diff > 1e-7 & nrep < 10) {
		
		nrep = nrep + 1
		evec = y - muhat    ##residuals a+e
		ecl = f_ecl(evec, ncl, szs)  ## residuals by group
        mod.lmer = NULL
		###estimate variance of a by lmer package### 
	#	mod.lmer = lmer(evec~-1+(1|id), REML=reml)
    #    thhat = summary(mod.lmer)$optinfo$val^2
		###estimate variance by grid search###
		ansi = try(ansi0<-uniroot(fth2rm, c(1e-10, 1e+3), szs=szs, ycl=ecl, N=n, xcl=xms, p=edf, type='ub', xtx=xtx, xtx2=xtx2, xmat_face=dd, ones=ones), silent=TRUE)
        if (class(ansi) == "try-error") {
            thhat = 0
        } else {
            thhat = ansi$root
        }
		
		type = "ub"
############update mu gaven a ############		
	ytil = NULL 
#gtil is edges
			gtil = NULL
			st = 1
			ed = 0
			sz = max(szs)
            pos = which(szs == sz)[1]
            oneMat = ones[[pos]]
			vi = diag(sz) + oneMat*thhat  ##covariance matrix
            covi = vi
            umat = t(chol(covi))
            uinv = solve(umat)
            #uinv0 is used for unbalanced
            uinv0 = uinv
			######L^-1*y=L^-1*(mu+xb)+e  e~N(0,I) #########
			for (icl in 1:ncl) {
				sz = szs[icl]
                uinv = uinv0[1:sz, 1:sz, drop=FALSE]
				yi = ycl[[icl]]
				ytil = c(ytil, uinv %*% as.matrix(yi, ncol=1))
				ed = ed + sz
				gtil = rbind(gtil, uinv %*% gmat[st:ed, ,drop=F])
				st = ed + 1
			}
			#####weighted coneB #########
			dsend = gtil[, (np + 1):(np + capm), drop = FALSE]
            zsend = gtil[, 1:np, drop = FALSE]
            ans = coneB(ytil, dsend, vmat = zsend, face=face)
            edf = ans$df
            face = ans$face
            bh = coef(ans)
			    if (any(round(bh[1:np],6) < 0)) {
                pos = (1:np)[which(round(bh[1:np],6) < 0)]
                face = unique(c(pos, face))
                }
			muhat = t(bigmat) %*% bh
			diff = mean((oldmu - muhat)^2)
			oldmu = muhat
            dd = t(bigmat[face, ,drop = FALSE])
            dd2 = gtil[,face,drop=FALSE]
                xms = list()
                st = 1
                ed = 0
                for (icl in 1:ncl) {
                    sz = szs[icl]
                    ed = ed + sz
                    xms[[icl]] = dd[st:ed, ,drop=F]
                    st = ed + 1
                }
			}
		ebars = sapply(ecl, mean)
		sig2hat = fsig(thhat, szs, ecl, ncl, N=n, edf=edf, D=nrow(bigmat), type=type)
		siga2hat = sig2hat * thhat 
		ahat = ebars*szs*thhat/(1+szs*thhat)
		rslt = list(muhat = muhat,  bh = bh,ahat = ahat, sig2hat = sig2hat, siga2hat = siga2hat, thhat = thhat, bigmat = bigmat,df=edf)
        return (rslt)
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

makedelta = function(x, sh, numknots = 0, knots = 0, space = "E", suppre = FALSE, interp = FALSE) {
#if (!interp) {
#x = (x - min(x)) / (max(x) - min(x))
#}
	n = length(x)
# find unique x values
#round(x,8) will make 0 edge in amat!
	#xu = sort(unique(round(x, 8)))
#new: center and scale to avoid numerical instabillity
	xu = sort(unique(x))
	n1 = length(xu)
	sm = 1e-7
	ms = NULL
#  increasing or decreasing
	if (sh < 3) {
		amat = matrix(0, nrow = n1 - 1, ncol = n)
		for (i in 1: (n1 - 1)) {
			amat[i, x > xu[i]] = 1
		}
		if (sh == 2) {amat = -amat}
		if (!interp) {
			for (i in 1:(n1 - 1)) {
#new: use ms in predict.cgam
				ms = c(ms, mean(amat[i, ]))
				amat[i, ] = amat[i, ] - mean(amat[i, ])
			}
		}
	} else if (sh == 3 | sh == 4) {
#  convex or concave
		amat = matrix(0, nrow = n1 - 2 ,ncol = n)
		#for (i in 1: (n1 - 2)) {
		#	amat[i, x > xu[i]] = x[x > xu[i]] - xu[i]
		#}
		for (i in 1: (n1 - 2)) {
			amat[i, x > xu[i+1]] = x[x > xu[i+1]] - xu[i+1]
		}
		if (sh == 4) {amat = -amat}
		xm = cbind(1:n*0+1,x)
		xpx = solve(t(xm) %*% xm)
		pm = xm %*% xpx %*% t(xm)
#new: use ms in predict.cgam
		if (!interp) {
			ms = amat %*% t(pm)
			#amat = amat - amat %*% t(pm)
			amat = amat - ms
		}
	} else if (sh > 4 & sh < 9) {
		amat = matrix(0, nrow = n1 - 1, ncol = n)
		if (sh == 5) { ### increasing convex
			for (i in 1:(n1 - 1)) {
				amat[i, x > xu[i]] = (x[x > xu[i]] - xu[i]) / (max(x) - xu[i])
			}
			if (!interp) {
				for (i in 1:(n1 - 1)) {
					ms = c(ms, mean(amat[i, ]))
					amat[i,] = amat[i,] - mean(amat[i,])
				}
			}
		} else if (sh == 6) {  ## decreasing convex
			for (i in 1:(n1 - 1)) {
				amat[i, x < xu[i + 1]] = (x[x < xu[i + 1]] - xu[i + 1]) / (min(x) - xu[i + 1])
			}
			if (!interp) {
				for (i in 1:(n1 - 1)) {
					ms = c(ms, mean(amat[i, ]))			
					amat[i,] = amat[i,] - mean(amat[i, ])
				}
			}
#print (ms)
		} else if (sh == 7) { ## increasing concave
			for (i in 1:(n1 - 1)) {
				amat[i, x < xu[i + 1]] = (x[x < xu[i + 1]] - xu[i + 1]) / (min(x) - xu[i + 1])
			}
			if (!interp) {
				for (i in 1:(n1 - 1)) {
					ms = c(ms, mean(amat[i, ]))
					amat[i,] = -amat[i,] + mean(amat[i,])
				}		
			}
		} else if (sh == 8) {## decreasing concave
			for (i in 1:(n1 - 1)) {
				amat[i, x > xu[i]] = (x[x > xu[i]] - xu[i]) / (max(x) - xu[i])
			}
			if (!interp) {
				for (i in 1:(n1 - 1)) {
					ms = c(ms, mean(amat[i, ]))
					amat[i,] = -amat[i,] + mean(amat[i,])
				}
			}
		}
	} else if (sh > 8 & sh < 18) {
        #new: add two knots
		#if (all(knots == 0) & numknots == 0) {
		if (length(knots) < 2 & numknots == 0) {
			if (sh == 9 | sh == 10) {#1 2
                #k = trunc(n1^(1/5)) + 4
                if (n1 <= 50) {
                    k = 5
                } else if (n1>50 && n1<100) {
                    k = 6
                } else if (n1>= 100 && n1<200) {
                    k = 7
                } else {
                    k = trunc(n1^(1/5)) + 6
                }
			} else {
                #k = trunc(n1^(1/7) + 4)
                if (n1 <= 50) {
                    k = 5
                } else if (n1>50 && n1<100) {
                    k = 6
                } else if (n1>= 100 && n1<200) {
                    k = 7
                } else {
                    k = trunc(n1^(1/7)) + 6
                }
            }
			if (space == "Q") {
				t = quantile(xu, probs = seq(0, 1, length = k), names = FALSE)
			}
			if (space == "E") {
				#t = 0:k / k * (max(x) - min(x)) + min(x)
				t = 0:(k-1) / (k-1) * (max(x) - min(x)) + min(x)
			} 
		#} else if (any(knots != 0) & numknots == 0) {
		} else if (length(knots) >= 2 & numknots == 0) {
			t = knots
		#} else if (all(knots == 0) & numknots != 0) {
		} else if (length(knots) < 2 & numknots != 0) {
			if (space == "Q") {
				t = quantile(xu, probs = seq(0, 1, length = numknots), names = FALSE)
			} 
			if (space == "E") {
				#k = numknots
#new: numknots should be the # of all knots
				k = numknots - 1
				#if (sh == 9 | sh == 10) {#1 2
				#	k = trunc(n1^(1/5)) + 4
				#} else {k = trunc(n1^(1/7) + 4)}
				t = 0:k / k * (max(x) - min(x)) + min(x)
			}
		#} else if (any(knots != 0) & numknots != 0) {
		} else if (length(knots) >= 2 & numknots != 0) {
			#t0 = quantile(xu, probs = seq(0, 1, length = numknots), names = FALSE)
			t = knots
			if (!suppre) {
				print("'knots' is used! 'numknots' is not used!")
			}
			#print ("'knots' is used!")
			#if (numknots != length(knots)) {
			#	if (!suppre) {
			#		print("length(knots) is not equal to 'numknots'! 'knots' is used!")
			#	}
			#} else if (any(t0 != knots)) {
			#	if (!suppre) {
			#		print("equal x-quantiles knots != 'knots'! 'knots' is used! ") 
			#	}
			#}
		}
		if (sh == 9) {#1			
			amat_ans = monincr(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 10) {#2
			amat_ans = mondecr(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 11) {#3
			amat_ans = convex(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 12) {#4
			amat_ans = concave(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 13) {#5
			amat_ans = incconvex(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 14) {#6
			amat_ans = incconcave(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 15) {#7
			#amat_ans = -incconcave(x, t, interp)
			amat_ans = incconcave(x, t, interp)
			amat = -amat_ans$sigma
			if (!interp) {
				ms = -amat_ans$ms
			}
		} else if (sh == 16) {#8
			#amat_ans = -incconvex(x, t, interp)
			amat_ans = incconvex(x, t, interp)
			amat = -amat_ans$sigma
			if (!interp) {
				ms = -amat_ans$ms
			}
		} else if (sh == 17) {#unconstrained
			amat_ans = incconvex(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
			#amat = -incconcave(x, t)
			#amat = rbind(x, t(bcspl(x, m = length(t), knots = t)$bmat)) 
			#amat = rbind(x, convex(x, t))
		}
	}
	#if (sh < 9) {
	#	rslt = list(amat = amat, knots = 0, ms = ms)
	#} else {
	#	rslt = list(amat = amat, knots = t, ms = ms)
	#}
	if (sh < 9) {t = 0}	
	rslt = list(amat = amat, knots = t, ms = ms)
	rslt
}

# Monotone increasing
monincr = function(xs, t, interp = FALSE) {
	n = length(xs)
#xs = (xs - min(xs)) / (max(xs) - min(xs))
	x = sort(xs)
	k = length(t) - 2
	m = k + 2
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	knt = 1:m
	for (i in 1:(k+2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	for (j in 1:(k-1)) {
		index = x >= t[1] & x <= t[j]
		sigma[j, index] = 0

		index = x > t[j] & x <= t[j+1]
		sigma[j, index] = (x[index] - t[j])^2 / (t[j+2] - t[j]) / (t[j+1] - t[j])

		index = x > t[j+1] & x <= t[j+2]
		sigma[j, index] = 1 - (x[index] - t[j+2])^2 / (t[j+2] - t[j+1]) / (t[j+2] - t[j])
	    
		index = x > t[j+2] #& x <= t[m]
		sigma[j, index] = 1
	}
	index = x >= t[1] & x <= t[k]
	sigma[k, index] = 0
	
	index = x > t[k] & x <= t[k+1]
	sigma[k, index] = (x[index] - t[k])^2 / (t[k+2] - t[k]) / (t[k+1] - t[k])
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k, index] = 1 - (x[index] - t[k+2])^2 / (t[k+2] - t[k+1]) / (t[k+2] - t[k])
	
	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = 1 - (t[2] - x[index])^2 / (t[2] - t[1])^2

	index = x > t[2] 
	sigma[k+1, index] = 1
	
	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = 0
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = (x[index] - t[k+1])^2 / (t[k+2] - t[k+1])^2
	
#new:
	ms = NULL
	if (!interp) {
		ms = apply(sigma, 1, mean)
		for (i in 1:m) {
			sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	}
	rslt = list(sigma = sigma, ms = ms)
	rslt
}	


########################################################
# Monotone decreasing
mondecr = function(xs, t, interp = FALSE) {
#xs = (xs - min(xs)) / (max(xs) - min(xs))
	x = sort(xs)
	n = length(x)
	k = length(t) - 2
	m = k + 2
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	#knt = 1:m
	#for (i in 1:(k + 2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	#t = x[knt]
	for (j in 1:(k - 1)) {
	 	index = x >= t[1] & x <= t[j]
	 	sigma[j, index] = 1

		index = x > t[j] & x <= t[j+1]
	 	sigma[j, index] = 1 - (x[index] - t[j])^2 / (t[j+2] - t[j]) / (t[j+1] - t[j])

	    	index = x > t[j+1] & x <= t[j+2]
	    	sigma[j, index] = (x[index] - t[j+2])^2 / (t[j+2] - t[j+1]) / (t[j+2] - t[j])

	    	index = x > t[j+2] 
	    	sigma[j, index] = 0
	}

	index = x >= t[1] & x <= t[k]
	sigma[k, index] = 1
	
	index = x > t[k] & x <= t[k+1]
	sigma[k, index] = 1 - (x[index] - t[k])^2 / (t[k+2] - t[k]) / (t[k+1] - t[k])

	index = x > t[k+1] & x <= t[k+2]
	sigma[k, index] = (x[index] - t[k+2])^2 / (t[k+2] - t[k+1]) / (t[k+2] - t[k])

	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = (t[2] - x[index])^2 / (t[2] - t[1])^2

	index = x > t[2] 
	sigma[k+1, index] = 0

	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = 1
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = 1 - (x[index] - t[k+1])^2 / (t[k+2] - t[k+1])^2

	ms = NULL
	if (!interp) {
		ms = apply(sigma, 1, mean)
		for (i in 1:m) {
			sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	}
	rslt = list(sigma = sigma, ms = ms)
	rslt
}

########################################################
# Convex
convex = function(xs, t, interp = FALSE) {
#xs = (xs - min(xs)) / (max(xs) - min(xs))
	x = sort(xs)
	n = length(x)
	k = length(t) - 2
	m = k + 2
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	#knt = 1:m
	#for (i in 1:(k+2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	for (j in 1:(k-1)) {
	 	index = x >= t[1] & x <= t[j]
	 	sigma[j, index] = 0
	 	
	 	index = x > t[j] & x <= t[j+1]
	 	sigma[j, index] = (x[index] - t[j])^3 / (t[j+2] - t[j]) / (t[j+1] - t[j]) / 3
	    
	   	index = x > t[j+1] & x <= t[j+2]
	    	sigma[j, index] = x[index] - t[j+1] - (x[index] - t[j+2])^3 / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) / 3 + (t[j+1] - t[j])^2 / 3 /(t[j+2] - t[j]) - (t[j+2] - t[j+1])^2 / 3 / (t[j+2] - t[j])

	    	index = x > t[j+2]
	    	sigma[j, index] = (x[index] - t[j+1]) + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) - (t[j+2] - t[j+1])^2 / 3 / (t[j+2] - t[j])
	}
	index = x >= t[1] & x <= t[k]
	sigma[k, index] = 0
	
	index = x > t[k] & x <= t[k+1]
	sigma[k, index] = (x[index] - t[k])^3 / (t[k+2] - t[k]) / (t[k+1] - t[k]) / 3

	index = x > t[k+1] & x <= t[k+2]
	sigma[k, index] = x[index] - t[k+1] - (x[index] - t[k+2])^3 / (t[k+2] - t[k]) / (t[k+2] - t[k+1]) / 3 + (t[k+1] - t[k])^2 / 3 / (t[k+2] -t[k]) - (t[k+2] - t[k+1])^2 / 3 / (t[k+2] - t[k])
	
	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = x[index] - t[1] + (t[2] - x[index])^3 / (t[2] - t[1])^2 / 3 - (t[2] - t[1]) / 3 #-(t[2]-t[1])^3/(t[2]-t[1])^2/3 #
	
	index = x > t[2] 
	sigma[k+1, index] = x[index] - t[1] - (t[2] - t[1]) / 3 #-(t[2]-t[1])^3/(t[2]-t[1])^2/3#
		
	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = 0
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = (x[index] - t[k+1])^3 / (t[k+2] - t[k+1])^2 / 3
	
	ms = NULL
	if (!interp) {
		xm = cbind(1:n*0+1, x)
		pm = xm %*% solve(t(xm) %*% xm) %*% t(xm)
		ms = matrix(0, nrow = nrow(sigma), ncol = ncol(sigma))
		for (i in 1:m) {
			ms[i,] = pm %*% sigma[i,]
			ms[i,] = ms[i, rank(xs)]		
			#rng=max(sigma[i,])
			#sigma[i,]=sigma[i,]/rng
			sigma[i,] = sigma[i,] - pm %*% sigma[i,]
			sigma[i,] = sigma[i, rank(xs)]
		}
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - pm %*% sigma[i,]
			sigma[i,] = sigma[i, rank(xs)]
		}
	}
	rslt = list(sigma = sigma, ms = ms)
	rslt
}


########################################################
# Concave
concave = function(xs, t, interp = FALSE) {
#xs = (xs - min(xs)) / (max(xs) - min(xs))
	x = sort(xs)
	n = length(x)
	k = length(t) - 2
	m = k + 2
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	#knt = 1:m
	#for (i in 1:(k+2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	#t = x[knt]
	for (j in 1:k) {
	 	index = x >= t[1] & x <= t[j]
	 	sigma[j, index] = x[index] - t[1]
	 	
	 	index = x > t[j] & x <= t[j+1]
	 	sigma[j, index] = t[j] - t[1] + ((t[j+1] - t[j])^3 - (t[j+1] - x[index])^3) / 3 / (t[j+1] - t[j]) / (t[j+2] - t[j]) + (x[index] - t[j]) * (t[j+2] - t[j+1]) / (t[j+2] - t[j])
	    
	        index = x > t[j+1] & x <= t[j+2]
	    	sigma[j, index] = t[j] - t[1] + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) + (t[j+2] - t[j+1]) * (t[j+1] - t[j]) / (t[j+2] - t[j]) + ((t[j+2] - t[j+1])^3 - (t[j+2] - x[index])^3) / 3 / (t[j+2] - t[j+1]) / (t[j+2] - t[j])	
 	   
 	   	index = x > t[j+2]
 	   	sigma[j, index] = t[j] - t[1] + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) + (t[j+2] - t[j+1]) * (t[j+1] - t[j]) / (t[j+2] - t[j]) + (t[j+2] - t[j+1])^2 / 3 / (t[j+2] - t[j])
	}

	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = -(t[2] - x[index])^3 / 3 / (t[2] - t[1])^2
	
	index = x > t[2] 
	sigma[k+1, index] = 0
	
	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = x[index] - t[1]
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = t[k+1] - t[1] + ((t[k+2] - t[k+1])^2 * (x[index] - t[k+1]) - (x[index] - t[k+1])^3 / 3) / (t[k+2] - t[k+1])^2
	
	ms = NULL
	if (!interp) {
		xm = cbind(1:n*0+1, x)
		pm = xm %*% solve(t(xm) %*% xm) %*% t(xm)
		ms = matrix(0, nrow = nrow(sigma), ncol = ncol(sigma))
		for (i in 1:m) {
			ms[i,] = pm %*% sigma[i,]
			ms[i,] = ms[i, rank(xs)]		
			sigma[i,] = sigma[i,] - pm %*% sigma[i,]
			sigma[i,] = sigma[i, rank(xs)]
		}
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - pm %*% sigma[i,]
			sigma[i,] = sigma[i, rank(xs)]
		}
	}	
	rslt = list(sigma = sigma, ms = ms)
	rslt
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
#########################################
#             CODE END                  #
#########################################