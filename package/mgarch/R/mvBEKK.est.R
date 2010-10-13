## Copyright (C) 2004 Harald SCHMIDBAUER - Vehbi Sinan TUNALIOGLU
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.


####################################################
## TODO LIST
## 1. t-distribution is to be added 
####################################################

mvBEKK.est<-
function(
		eps,			# an data frame holding the time series
		order  = c(1,1),	# order of the BEKK(p,q) model to be estimated c(p,q)
		params = NULL,		# initial parameters for the optim function
		fixed  = NULL,		# parameter list that are to be fixed
		method = "BFGS",	# the method that will be used in the optim process
		verbose = F
	)
{
  
  count.triangular <-
    function(dimension){
      if(dimension <= 0){
        0
      }
      else{
        dimension + count.triangular(dimension - 1)
      }
    }
  
  if(verbose == T){
    out <- function(...){
      cat(...)
    }
  }
  else{
    out <- function(...) { }
  }
  
    # get the length and the number of the series
	series.length = length(eps[,1])
	series.count  = length(eps[1,])

	# check the given order
	# orders should be integers
	if(order[1] != as.integer(order[1]) || order[2] != as.integer(order[2]))
	{
		stop("order property should contain integer values")
	}
	
	# GARCH effect could be set to 0, but, ARCH should be greater than 0
	if(order[1] < 0 || order[2] < 1)
	{
		stop("BEKK(",order[1],",",order[2],") is not implemented.")
	}
	
	# construct the paramters list.
	# first get the length of the parameter list
	params.length = count.triangular(series.count) + (order[2] * series.count^2) + (order[1] * series.count^2) 

	if(is.null(params))
	{
		# TODO
		# these are meaningless parameters.
		# set some useful initial parameters.
		params = c(1,0,1,0,0,1)
		params = c(params, rep(0.1, params.length - 6))
		out("\nWarning: initial values for the parameters are set to:\n\t", params,"\n")
	}
	else if(length(params) != params.length)
	{
		stop("Length of the initial parameter list doesn't conform required length. There should be ", params, " parameters in total")
	}

	# check the given fixed parameters
	if(!is.null(fixed))
	{
		# check the format of the fixed parameters
		if(
				(!is.array(fixed)) ||
				(dim(fixed)[1] != 2) ||
				(length(fixed[1,]) != length(fixed[2,])))
		{
			stop("fixed should be an array of two vectors. Try fixed = array(c(a,b,c,d,...), dim = c(2,y))")
		}

		# check the first dimension, if it contains appropriate index values,
		# that is integer values rather than floating or negative numbers
		for(count in 1:length(fixed[1,]))
		{
			if((fixed[1,count] != as.integer(fixed[1,count])) || (fixed[1,count] <= 0))
			{
				stop("First dimension of the fixed array should contain only positive integer values for indexing purposes")
			}
		}

		# check the length of the fixed parameters
		if(length(fixed[1,]) > length(params))
		{
			stop("fixed array could not contain more index-value pairs than the params array length");
		}
	}
	
	# check the method specified in the argument list
	if(!(
		(method == "Nelder-Mead") ||
		(method == "BFGS") ||
		(method == "CG") ||
		(method == "L-BFGS-B") ||
		(method == "SANN")
		))
	{
		stop("'", method, "' method is not available")
	}

	
	fake.params = params
	if(!is.null(fixed))
	{
		# extract the parameters specified in the fixed list.
		fake.params = params
		for(i in 1:length(fixed[1,]))
		{
			fake.params[fixed[1,][i]] = NA
		}
		fake.params = na.omit(fake.params)
	}
	
	# parameters seem appropriate
	# define the loglikelihood function
	loglikelihood.C <- function(params)
	{
		loglikelihood.C <- .C("loglikelihood",
			as.vector(params, mode = "double"),
			as.vector(fixed[1,], mode = "integer"),
			as.vector(fixed[2,], mode = "double"),
			as.integer(length(fixed[1,])),
			as.vector(t(eps)), # funny: transpose the time series
			as.integer(series.count),
			as.integer(series.length),
			as.vector(order, mode = "integer"),
			retval = 0.0,
			PACKAGE = "mgarch"
		)
		
		if(is.nan(loglikelihood.C$retval) == T)
		{
			nonusedret = 1e+100
			
		}
		else
		{
			nonusedret = loglikelihood.C$retval
		}
		nonusedret
	}

	# begin estimation process
	
	# first log the start time
	start = Sys.time()		
	out("* Starting estimation process.\n")
	out("* Optimization Method: '", method, "'\n")

	# call the optim function
	estimation = optim(fake.params, loglikelihood.C, method = method, hessian = T)

	# estimation completed
	out("* Estimation process completed.\n")

	# log estimation time
	est.time = difftime(Sys.time(), start)

	# calculate the AIC
	# it is estimation value + number of estimated parameters (punishment :))
	aic = estimation$value + (length(params) - length(fixed[1,]))

	# following script will prepare an object that holds the estimated
	# parameters and some useful diagnostics data like estimated correlation,
	# standard deviation, eigenvalues and so on.

	# TODO
	# estimation$hessian is non-existing if fixed parameter list contains all the
	# paramters to be estimated. That is that the estimation procedure gets no parameters,
	# thus, there is no errors... Fix it... How? 
	# Whether encapsulate with an "if" statement, probably not efficient,
	# or give a fake hessian

	# give a fake hessian
	if(length(fake.params) == 0)
	{
		estimation$hessian = matrix(rep(0.1, series.count^2), nrow = series.count, ncol = series.count)
	}

	# get the hessian matrix and grap the diagonal
	inv.hessian.mat = solve(estimation$hessian)

	diag.inv.hessian = sqrt(abs(diag(inv.hessian.mat)))
	if(length(which(diag(inv.hessian.mat) < 0)) == 0)
	{
		warning("negative inverted hessian matrix element")
	}
	
	# fix the asymptotic-theory standard errors of the
	# coefficient estimates with fixed parameters
	if(!is.null(fixed))
	{
		temp.diag.inv.hessian = numeric()
		shifted = 0
		for(count in 1:params.length)
		{
			check.point = 0
			for(i in 1:length(fixed[1,]))
			{
				if(count == fixed[1,i])
				{
					check.point = 1
					shifted = shifted + 1
					temp.diag.inv.hessian[count] = 0
					break
				}
			}
			if(check.point == 0)
			{
				temp.diag.inv.hessian[count] = diag.inv.hessian[count - shifted]
			}
		}
		diag.inv.hessian = temp.diag.inv.hessian
	}
	# construct the asymptotic-theory standard errors of the coefficient estimates matrices
	parnum = 1 + order[1] + order[2]	# calculate number of paramater matrices
	asy.se.coef = list()			# declare the asy.se.coef matrices list

	############################################################################
	### CRITICAL!!!
	### since we are not bidimensional anymore, be careful!!!
	############################################################################
	# first initialize the first asy.se.coef matrix, corresponding to the C matrix
	tmp.array = array(rep(0, series.count^2), dim = c(series.count, series.count))
	tmp.array[!lower.tri(tmp.array)] = diag.inv.hessian[1:length(which(!lower.tri(tmp.array) == T))]
	asy.se.coef[[1]] = tmp.array 

	# following loop initalizes the ARCH and GARCH parameter matrices respectively
	for(count in 1:(parnum - 1))
	{
		# !! a bit hard to follow
		asy.se.coef[[count + 1]] = array(diag.inv.hessian[(count.triangular(series.count) + 1 + (count - 1) * series.count^2):(count.triangular(series.count) + 1 + series.count^2 + (count - 1) * series.count^2)], dim = c(series.count, series.count));
	}
	
	buff.par = list()		# declare the parameter list
	
	# shift the fixed parameters inside the estimated parameters
	if(!is.null(fixed))
	{
		estim.params = numeric()
		shifted = 0
		for(count in 1:params.length)
		{
			check.point = 0
			for(i in 1:length(fixed[1,]))
			{
				if(count == fixed[1,i])
				{
					check.point = 1
					shifted = shifted + 1
					estim.params[count] = fixed[2,i]
					break
				}
			}
			if(check.point == 0)
			{
				estim.params[count] = estimation$par[count - shifted]
			}
		}
	}
	else
	{
		estim.params = estimation$par
	}

	# first initialize the C matrix
	tmp.array = array(rep(0, series.count^2), dim = c(series.count, series.count))
	tmp.array[!lower.tri(tmp.array)] = estim.params[1:length(which(!lower.tri(tmp.array) == T))]

	buff.par[[1]] = tmp.array

	# following loop initalizes the ARCH and GARCH parameter matrices respectively
	for(count in 1:(parnum - 1))
	{
		# !! a bit hard to follow
		buff.par[[count + 1]] = array(estim.params[(count.triangular(series.count) + 1 + (count - 1) * series.count^2):(count.triangular(series.count) + 1 + series.count^2 + (count - 1) * series.count^2)], dim = c(series.count, series.count));
	}

	# calculate the transposes of the parameter matrices
	buff.par.transposed = lapply(buff.par, t)
	
	# start diagnostics
	out("* Starting diagnostics...\n")
	out("* Calculating estimated:\n")
	out("*\t1. residuals,\n")
	out("*\t2. correlations,\n")
	out("*\t3. standard deviations,\n")
	out("*\t4. eigenvalues.\n")

	HLAGS = list()		# list of H lags that will be used later in the MGARCH implementation
	for(count in 1:order[1])
	{
		# TODO:check intial values (currently 1's on the diagonal)
		HLAGS[[count]] = array(rep(0, series.count^2), dim = c(series.count,series.count))
		diag(HLAGS[[count]]) = 1
	}

	residuals = list()
	for(i in 1:series.count)
	{
		residuals[[i]] = numeric()
	}

	# initialize the first residuals we are not able to calculate
	for(count in 1:max(order))
	{
		for(i in 1:series.count)
		{
			residuals[[i]][count] = 0
		}
	}
	
	resid = array(rep(0,series.count), dim = c(series.count,1)) # declare a temporary residuals buffer

	# calculate eigenvalues
	# TODO:
	# Angi says that following is not true according to Bauwens, Laurent, Rombouts Paper.
	temp = 0
	for(count in 2:parnum)
	{
		temp = temp + kronecker(buff.par[[count]], buff.par[[count]]) 
	}
	eigenvalues = svd(temp)$d

	##################################################################
	### TODO:
	### FROM NOW ON, HELP NEEDED
	### ASK HARALD HOCA!
	##################################################################

	# compute the unconditional covariance matrix
	numerat = t(buff.par[[1]]) %*% buff.par[[1]]
	dim(numerat) = c(series.count^2,1)
	denom = solve(diag(rep(1, series.count^2)) - temp)
	sigma = denom %*% numerat
	dim(sigma) = c(series.count, series.count)
	
	H = cov(eps)	# to initialize, use the covariance matrix of the series 
	H.estimated = lapply(1:series.length, function(x){H})

	cor = list()		# declare the estimated correlation series
	for(i in 1:series.count)
	{
		cor[[i]] = list()
		for(j in 1:series.count)
		{
			cor[[i]][[j]] = numeric()
		}
	}
	sd = list()		# declare the estimated standard deviation series
	for(i in 1:series.count)
	{
		sd[[i]] = numeric()
	}
	
	eps.est = array(rep(0,series.count), dim = c(series.count,1))	# declare a temporary eps buffer
	
	CTERM = buff.par.transposed[[1]] %*% buff.par[[1]] # calculate the C'C term

	out("* Entering Loop...");
	for(count in (max(order) + 1):series.length) # cruical loop! initializing diagnostics data
	{
		# do the swap calculation for H terms
		if(order[1] >= 2)
		{
			for(tmp.count in order[1]:2)
			{
				HLAGS[[tmp.count]] = HLAGS[[(tmp.count - 1)]]
			}
		}
		HLAGS[[1]] = H

		# a bit complicated but following explanation will be useful hopefully
		# H = (C')x(C) + (A')(E_t-1)(E_t-1')(A) + (B')(E_t-2)(E_t-2')(B) + ... +  (G')(H_t-1)(G) + (L')(H_t-2)(L) + ... 
		#                    |_____________|          |_____________|             |____________|   |____________| |_____|
		#                        E1 TERM                  E2 TERM                     G1 TERM         G2 TERM     G3.G4..
		#                |____________________|   |____________________| |_____|
		#                        A1 TERM                  A2 TERM        A3.A4..
		#	   |______|  |_____________________________________________________|  |______________________________________|  
		#      C TERM                         A TERM                                              G TERM
		
		H = CTERM
		ord1 = 1
		for(tmp.count in 1:(order[2] + order[1]))
		{
			if(tmp.count <= order[2])
			{
				# ARCH EFFECT (A TERM)
				H = H + buff.par.transposed[[tmp.count + 1]] %*% as.matrix(t(eps[count - tmp.count,])) %*% as.matrix(eps[count - tmp.count,]) %*% buff.par[[tmp.count + 1]]
			}
			else
			{
				# GARCH EFFECT (G TERM)
				H = H + buff.par.transposed[[tmp.count + 1]] %*% HLAGS[[ord1]]  %*% buff.par[[tmp.count + 1]]
				ord1 = ord1 + 1
			}
		}

		# TODO add appropriate comments for following assignments and calculations
  		H.estimated[[count]] = H 
		svdH = svd(H)
		sqrtH = svdH$u %*% diag(sqrt(svdH$d)) %*% t(svdH$v)
		
		invsqrtH = solve(sqrtH)
		resid = invsqrtH %*% as.matrix(t(eps[count,]))
		for(i in 1:series.count)
		{
			residuals[[i]][count] = resid[i,1]
		}
	
		# TODO: check
		for(i in 1:series.count)
		{
			for(j in 1:series.count)
			{
				cor[[i]][[j]][count] = H[i,j] / sqrt(H[i,i] * H[j,j])
			}
		}
		for(i in 1:series.count)
		{
			sd[[i]][count] = sqrt(H[i,i])
		}
	}
	
	# diagnostics ready
	out("Diagnostics ended...\n")
	
	names(order) <- c("GARCH component", "ARCH component")
	names(buff.par) <- as.integer(seq(1, parnum))

	mvBEKK.est <- list(
		eps = eps,
		series.length = series.length,
		estimation.time = est.time,
		total.time = difftime(Sys.time(), start),
		order = order,
		estimation = estimation,
		aic = aic,
		asy.se.coef = asy.se.coef,
		est.params = buff.par,
		cor = cor,
		sd = sd,
		H.estimated = H.estimated,
		eigenvalues = eigenvalues,
		uncond.cov.matrix = sigma,
		residuals = residuals
	)

	class(mvBEKK.est) = "mvBEKK.est"
	
	out("Class attributes are accessible through following names:\n")
	out(names(mvBEKK.est), "\n")
	
	return(mvBEKK.est)
}
