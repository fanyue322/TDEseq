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

convex.lm <-function(x,t){
	n=length(x)
	k=length(t)-2
	m=k+2
	sigma=matrix(1:m*n,nrow=m,ncol=n)
	dsigma=matrix(1:m*n,nrow=m,ncol=n)
	for(j in 1:(k-1)){
		i1=x<=t[j]
		sigma[j,i1] = 0
		dsigma[j,i1] = 0
	 	i2=x>t[j]&x<=t[j+1]
	 	sigma[j,i2] = (x[i2]-t[j])^3 / (t[j+2]-t[j]) / (t[j+1]-t[j])/3
	 	dsigma[j,i2] = 3*(x[i2]-t[j])^2 / (t[j+2]-t[j]) / (t[j+1]-t[j])/3
	    i3=x>t[j+1]&x<=t[j+2]
	    sigma[j,i3] = x[i3]-t[j+1]-(x[i3]-t[j+2])^3/(t[j+2]-t[j])/(t[j+2]-t[j+1])/3+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
	    dsigma[j,i3] = 1-3*(x[i3]-t[j+2])^2/(t[j+2]-t[j])/(t[j+2]-t[j+1])/3
	    i4=x>t[j+2]
	    sigma[j,i4]=(x[i4]-t[j+1])+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
	    dsigma[j,i4]=1
	}
	i1=x<=t[k]
	sigma[k,i1] = 0
	dsigma[k,i1] = 0
	i2=x>t[k]&x<=t[k+1]
	sigma[k,i2] = (x[i2]-t[k])^3 / (t[k+2]-t[k]) / (t[k+1]-t[k])/3
	dsigma[k,i2] = 3*(x[i2]-t[k])^2 / (t[k+2]-t[k]) / (t[k+1]-t[k])/3
	i3=x>t[k+1]
	sigma[k,i3] = x[i3]-t[k+1]-(x[i3]-t[k+2])^3/(t[k+2]-t[k])/(t[k+2]-t[k+1])/3+(t[k+1]-t[k])^2/3/(t[k+2]-t[k])-(t[k+2]-t[k+1])^2/3/(t[k+2]-t[k])
	dsigma[k,i3] = 1-3*(x[i3]-t[k+2])^2/(t[k+2]-t[k])/(t[k+2]-t[k+1])/3
	i1=x<=t[2]
	sigma[k+1,i1]=x[i1]-t[1]+(t[2]-x[i1])^3/(t[2]-t[1])^2/3-(t[2]-t[1])^3/(t[2]-t[1])^2/3
	dsigma[k+1,i1]=1-3*(t[2]-x[i1])^2/(t[2]-t[1])^2/3
	i2=x>t[2]
	sigma[k+1,i2]=x[i2]-t[1]-(t[2]-t[1])^3/(t[2]-t[1])^2/3
	dsigma[k+1,i2]=1
	i1=x<=t[k+1]
	sigma[k+2,i1]=0
	dsigma[k+2,i1]=0
	i2=x>t[k+1]
	sigma[k+2,i2]=(x[i2]-t[k+1])^3/(t[k+2]-t[k+1])^2/3
	dsigma[k+2,i2]=3*(x[i2]-t[k+1])^2/(t[k+2]-t[k+1])^2/3
	
	for(i in 1:m){
		rng=max(sigma[i,])
		sigma[i,]=sigma[i,]/rng
		dsigma[i,]=dsigma[i,]/rng
	}
	
	ans=new.env()
	ans$sigma=sigma
	ans$dsigma=dsigma
	ans
}

concave.lm <-function(x,t){
	n=length(x)
	k=length(t)-2
	m=k+2
	sigma=matrix(1:m*n,nrow=m,ncol=n)
	dsigma=matrix(1:m*n,nrow=m,ncol=n)
	for(j in 1:(k-1)){
		i1=x<=t[j]
		sigma[j,i1] = 0
		dsigma[j,i1] = 0
	 	i2=x>t[j]&x<=t[j+1]
	 	sigma[j,i2] = (x[i2]-t[j])^3 / (t[j+2]-t[j]) / (t[j+1]-t[j])/3
	 	dsigma[j,i2] = 3*(x[i2]-t[j])^2 / (t[j+2]-t[j]) / (t[j+1]-t[j])/3
	    i3=x>t[j+1]&x<=t[j+2]
	    sigma[j,i3] = x[i3]-t[j+1]-(x[i3]-t[j+2])^3/(t[j+2]-t[j])/(t[j+2]-t[j+1])/3+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
	    dsigma[j,i3] = 1-3*(x[i3]-t[j+2])^2/(t[j+2]-t[j])/(t[j+2]-t[j+1])/3
	    i4=x>t[j+2]
	    sigma[j,i4]=(x[i4]-t[j+1])+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
	    dsigma[j,i4]=1
	}
	i1=x<=t[k]
	sigma[k,i1] = 0
	dsigma[k,i1] = 0
	i2=x>t[k]&x<=t[k+1]
	sigma[k,i2] = (x[i2]-t[k])^3 / (t[k+2]-t[k]) / (t[k+1]-t[k])/3
	dsigma[k,i2] = 3*(x[i2]-t[k])^2 / (t[k+2]-t[k]) / (t[k+1]-t[k])/3
	i3=x>t[k+1]
	sigma[k,i3] = x[i3]-t[k+1]-(x[i3]-t[k+2])^3/(t[k+2]-t[k])/(t[k+2]-t[k+1])/3+(t[k+1]-t[k])^2/3/(t[k+2]-t[k])-(t[k+2]-t[k+1])^2/3/(t[k+2]-t[k])
	dsigma[k,i3] = 1-3*(x[i3]-t[k+2])^2/(t[k+2]-t[k])/(t[k+2]-t[k+1])/3
	i1=x<=t[2]
	sigma[k+1,i1]=x[i1]-t[1]+(t[2]-x[i1])^3/(t[2]-t[1])^2/3
	dsigma[k+1,i1]=1-3*(t[2]-x[i1])^2/(t[2]-t[1])^2/3
	i2=x>t[2]
	sigma[k+1,i2]=x[i2]-t[1]
	dsigma[k+1,i2]=1
	i1=x<=t[k+1]
	sigma[k+2,i1]=0
	dsigma[k+2,i1]=0
	i2=x>t[k+1]
	sigma[k+2,i2]=(x[i2]-t[k+1])^3/(t[k+2]-t[k+1])^2/3
	dsigma[k+2,i2]=3*(x[i2]-t[k+1])^2/(t[k+2]-t[k+1])^2/3
	
	for(i in 1:m){
		sigma[i,]=max(dsigma[i,])*x-sigma[i,]
		dsigma[i,]=max(dsigma[i,])-dsigma[i,]
	}
	for(i in 1:m){
		rng=max(sigma[i,])
		sigma[i,]=sigma[i,]/rng
		dsigma[i,]=dsigma[i,]/rng
	}
	
	
	ans=new.env()
	ans$sigma=sigma
	ans$dsigma=dsigma
	ans
}

monincr.lm <-function(x,t){
	n=length(x)
	k=length(t)-2
	m=k+2
	sigma=matrix(1:m*n,nrow=m,ncol=n)
	dsigma=sigma
	for(j in 1:(k-1)){
	 	i1=x<=t[j]
	 	sigma[j,i1] = 0
	 	dsigma[j,i1] = 0
	 	i2=x>t[j]&x<=t[j+1]
		sigma[j,i2] = (x[i2]-t[j])^2 / (t[j+2]-t[j]) / (t[j+1]-t[j])
		dsigma[j,i2] = 2*(x[i2]-t[j]) / (t[j+2]-t[j]) / (t[j+1]-t[j])
	    i3=x>t[j+1]&x<=t[j+2]
		sigma[j,i3] = 1-(x[i3]-t[j+2])^2/(t[j+2]-t[j+1])/(t[j+2]-t[j])
		dsigma[j,i3] = -2*(x[i3]-t[j+2])/(t[j+2]-t[j+1])/(t[j+2]-t[j])
	    i4=x>t[j+2]
		sigma[j,i4]=1
		dsigma[j,i4]=0
	}
	
	i1=x<=t[k]
	sigma[k,i1] = 0
	dsigma[k,i1] = 0
	i2=x>t[k]&x<=t[k+1]
	sigma[k,i2] = (x[i2]-t[k])^2 / (t[k+2]-t[k]) / (t[k+1]-t[k])
	dsigma[k,i2] = 2*(x[i2]-t[k]) / (t[k+2]-t[k]) / (t[k+1]-t[k])
	i3=x>t[k+1]&x<=t[k+2]
	sigma[k,i3] = 1- (x[i3]-t[k+2])^2/(t[k+2]-t[k+1])/(t[k+2]-t[k])
	dsigma[k,i3] = -2*(x[i3]-t[k+2])/(t[k+2]-t[k+1])/(t[k+2]-t[k])
	i4=x>t[k+2]
	sigma[k,i4]=1
	dsigma[k,i4]=0
	
	i1=x<=t[2]
	sigma[k+1,i1]=1-(t[2]-x[i1])^2/(t[2]-t[1])^2
	dsigma[k+1,i1]=2*(t[2]-x[i1])/(t[2]-t[1])^2
	i2=x>t[2]
	sigma[k+1,i2]=1
	dsigma[k+1,i2]=0
	
	i1=x<=t[k+1]
	sigma[k+2,i1]=0
	dsigma[k+2,i1]=0
	i2=x>t[k+1]&x<=t[k+2]
	sigma[k+2,i2]=(x[i2]-t[k+1])^2/(t[k+2]-t[k+1])^2
	dsigma[k+2,i2]=2*(x[i2]-t[k+1])/(t[k+2]-t[k+1])^2
	i3=x>t[k+2]
	sigma[k+2,i3]=1
	dsigma[k+2,i3]=0
	for(i in 1:m){
		rng=max(sigma[i,])
		sigma[i,]=sigma[i,]/rng
		dsigma[i,]=dsigma[i,]/rng
	}
	
	ans=new.env()
	ans$sigma=sigma
	ans$dsigma=dsigma
	ans
}