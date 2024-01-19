basisfunction<-function(x,knots=0,type,fit.model='lmm',nsim=100)
{
    if(fit.model!='lmm')
	{
	n=length(x)
    c=1.2
    one = 1:n * 0 + 1
    zmat = matrix(one, ncol = 1)
    k = 1
    add = 3
    if (type > 2) {
        add = 4
    }
    if (length(knots) > 1) {
        if (min(knots) <= min(x) & max(knots) >= max(x)) {
            t = knots
        }else {
            br = c(10, 25, 100, 200, 400, 1000, 1e+10)
            obs = 1:7
            nk = min(obs[n <= br]) + add
            t = 0:(nk - 1)/(nk - 1) * (max(x) - min(x)) + min(x)
        }
    }else {
        br = c(10, 25, 100, 200, 400, 1000, 1e+10)
        obs = 1:7
        nk = min(obs[n <= br]) + add
        t = 0:(nk - 1)/(nk - 1) * (max(x) - min(x)) + min(x)
    }
    if (type == 1) {
        bas = monincr.lm(x, t)
        delta = bas$sigma
        slopes = bas$dsigma
    }
    if (type == 2) {
        bas = monincr.lm(x, t)
        delta = 1 - bas$sigma
        slopes = -bas$dsigma
    }
    if (type == 3) {
        bas = convex.lm(x, t)
        delta = bas$sigma
        slopes = bas$dsigma
    }
    if (type == 4) {
        bas = concave.lm(x, t)
        delta = bas$sigma
        slopes = bas$dsigma
    }
    
	incr = 0
    decr = 0
    if (type == 1 | type == 5 | type == 7) {
        incr = 1
    }
    if (type == 2 | type == 6 | type == 8) {
        decr = 1
    }
    m = length(delta)/n
    if (incr == 0 & decr == 0) {
        zmat = cbind(zmat, x)
    }
	
	ztr = zmat
    dtr = delta
    
    k0 = dim(zmat)[2]
        mdist = 1:(m + 1) * 0
        for (isim in 1:nsim) {
            ysim = rnorm(n)
            asim = coneB(ysim, t(dtr), ztr)
            df0 = asim$df - k0
            mdist[df0 + 1] = mdist[df0 + 1] + 1
        }
        mdist = mdist/nsim
		
    return(list(delta=delta,slopes=slopes,mdist=mdist))
	}else{
	shape=type+8
	xmat=as.matrix(x)
    n = length(x)
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
        gmat = t(bigmat)
		
	    dsend = gmat[, (np + 1):(np + capm), drop = FALSE]
        zsend = gmat[, 1:np, drop = FALSE]
    
	
	
	    m=ncol(gmat)
		mdist=1:(m+1)*0
		k0=dim(zsend)[2]
		for(isim in 1:nsim){
			ysim=rnorm(n)
			asim=coneB(ysim,dsend,zsend)
			df0=asim$df-k0
			mdist[df0+1]=mdist[df0+1]+1
		}
		mdist=mdist/nsim
	}
	return(list(bigmat=bigmat,mdist=mdist))
}

conspline<-function(y,x,basis,type,zmat = 0, wt=0, test=TRUE, c=1.2, nsim=10000)
{
    delta=basis$delta
	mdist=basis$mdist
	slopes=basis$slopes
    n = length(y)
    if (n < 10) {
        print("ERROR: must have at least 10 observations")
    }
    if (length(wt) > 1 & min(wt) <= 0) {
        print("ERROR: must have positive weights")
    }
    one = 1:n * 0 + 1
    if (length(x) != length(y)) {
        print("ERROR: length of x must be length of y")
    }
    if (length(zmat) > 1) {
        if (length(zmat) == n) {
            zmat = matrix(zmat, ncol = 1)
        }
        if (dim(zmat)[1] != n) {
            print("ERROR: number of rows of zmat must be length of y")
        }
        k = dim(zmat)[2]
        rone = one - zmat %*% solve(t(zmat) %*% zmat) %*% t(zmat) %*% 
            one
        if (sum(rone^2) > 1e-08) {
            zmat = cbind(one, zmat)
            k = k + 1
        }
    }else {
        zmat = matrix(one, ncol = 1)
        k = 1
    }
    add = 3
    if (type > 2) {
        add = 4
    }
   
    incr = 0
    decr = 0
    if (type == 1 | type == 5 | type == 7) {
        incr = 1
    }
    if (type == 2 | type == 6 | type == 8) {
        decr = 1
    }
    m = length(delta)/n
    if (incr == 0 & decr == 0) {
        zmat = cbind(zmat, x)
    }
    if (length(wt) > 1) {
        ytr = y * sqrt(wt)
        dtr = delta
        ztr = zmat
        for (i in 1:n) {
            dtr[, i] = delta[, i] * sqrt(wt[i])
            ztr[i, ] = zmat[i, ] * sqrt(wt[i])
        }
    }else {
        ztr = zmat
        dtr = delta
        ytr = y
    }
    ans = coneB(ytr, t(dtr), ztr)
    dfuse = min(c * ans$df, m + k)
    sighat = sum((ytr - ans$yhat)^2)/(n - dfuse)
    if (k > 1) {
        use = abs(ans$coef) > 1e-08
        if (k > 1) {
            use[2:k] = FALSE
        }
        xj = cbind(ztr, t(dtr))
        xj = xj[, use]
        pj = xj %*% solve(t(xj) %*% xj) %*% t(xj)
        zm = ztr[, 2:k]
        ppinv = solve(t(zm) %*% (diag(one) - pj) %*% zm)
        zcoef = ans$coef[2:k]
        sez = sqrt(diag(ppinv) * sighat)
        tz = zcoef/sez
        pz = 2 * (1 - pt(abs(tz), n - dfuse))
    }
    if (length(wt) > 1) {
        muhat = ans$yhat/sqrt(wt)
    }
    else {
        muhat = ans$yhat
    }
    if (test) {
        th0 = ztr %*% solve(t(ztr) %*% ztr) %*% (t(ztr)%*%ytr)  ###will be more efficient for large scale data set
        sse0 = sum((ytr - th0)^2)
        sse1 = sum((ytr - ans$yhat)^2)
        bstat = (sse0 - sse1)/sse0
        k0=1
        ps = mdist[1]
        for (d in 1:m) {
            ps = ps + pbeta(bstat, d/2, (n - d - k0)/2) * mdist[d + 1]
        }
        pval = 1 - ps
	#	}else{
	#	d=ans$df
	#	ps = pbeta(bstat, d/2, (n - d - k0)/2)
	#	pval=1-ps
	#	}
    }
    if (incr == 0 & decr == 0) {
        fhat = t(delta) %*% ans$coef[(k + 2):(k + 1 + m)]
        fhat = fhat + x * ans$coef[k + 1]
        fhat = fhat + ans$coef[1]
        fslope = t(slopes) %*% ans$coef[(k + 2):(k + 1 + m)]
        fslope = fslope + ans$coef[k + 1]
    }
    else {
        fhat = t(delta) %*% ans$coef[(k + 1):(k + m)]
        fhat = fhat + ans$coef[1]
        fslope = t(slopes) %*% ans$coef[(k + 1):(k + m)]
    }
    cans = new.env()
    cans$sighat = sighat
    if (k > 1) {
        cans$zhmat = ppinv
        cans$zcoef = zcoef
        cans$sez = sqrt(sez)
        cans$pvalz = pz
    }
    if (test) {
        cans$pvalx = pval
    }
    cans$muhat = muhat
    cans$fhat = fhat
    cans$fslope = fslope
    cans$knots = t
    cans$df = dfuse
    wp1 = sum(ytr * ztr[, 1])/sum(ztr[, 1]^2) * ztr[, 1]
    cans$rsq = 1 - sum((ytr - ans$yhat)^2)/sum((ytr - wp1)^2)
    cans
}




############LMM function########################
conespline_lmm<-function(x,y,basis,group,shape=9,test=TRUE,nsim=100,mod.uniroot=mod.uniroot)
{ ## adjust
    bigmat=basis$bigmat
	mdist=basis$mdist
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

		if(mod.uniroot)
		{
		ansi = try(ansi0<-uniroot(fth2rm, c(1e-10, 1e+3), szs=szs, ycl=ecl, N=n, xcl=xms, p=edf, type='ub', xtx=xtx, xtx2=xtx2, xmat_face=dd, ones=ones), silent=TRUE)
        if (class(ansi) == "try-error") {
            thhat = 0
        } else {
            thhat = ansi$root
        }
		}else{
	#	 mod.lmer = try(mod.lmer0<-lmer(evec~-1+(1|id), REML=FALSE,verbose=0),silent=TRUE)
	    mod.lmer = try(mod.lmer0 <- nlme::lme(evec~ 1, random = ~1|id),silent=TRUE) 
		if(class(mod.lmer)=='try-error')
		{
		thhat=0
		}else{
		thhat = as.numeric(nlme::VarCorr(mod.lmer)[2,1])
		}
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
            th0 = zsend %*% solve(t(zsend) %*% zsend) %*% (t(zsend)%*%ytil)
			sse0=sum((ytil-th0)^2)
			sse1=sum((ytil-yhat)^2)
		bstat=(sse0-sse1)/sse0

		m=ncol(gtil)
		k0=dim(zsend)[2]
		# mdist=1:(m+1)*0
		# k0=dim(zsend)[2]
		# for(isim in 1:nsim){
			# ysim=rnorm(n)
			# asim=coneB(ysim,dsend,zsend)
			# df0=asim$df-k0
			# mdist[df0+1]=mdist[df0+1]+1
		# }
		# mdist=mdist/nsim
		ps=mdist[1]
		for(d in 1:m){
			ps=ps+pbeta(bstat,d/2,(n-d-k0)/2)*mdist[d+1]
		}
		pval=1-ps
    }
	rslt = list(muhat = muhat,  bh = bh,ahat = ahat, sig2hat = sig2hat, siga2hat = siga2hat, thhat = thhat, bigmat = bigmat,pval=pval,bstat=bstat)
    return (rslt)
}