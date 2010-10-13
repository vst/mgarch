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

mvBEKK.sim<- 
function(
			series.count,		# length of the series
			T,			# number of the series
			order  = c(1, 1),	# order list of the BEKK model : BEKK(p,q)
			params = NULL		# initial parameters for simulation 
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
  
	# check the given order
	# orders should be integers
	if(order[1] != as.integer(order[1]) || order[2] != as.integer(order[2]))
	{
		stop("order should contain integer values")
	}
	# GARCH effect could be set to 0, but, ARCH could not be 0
	if(order[1] < 0 || order[2] < 1)
	{
		stop("BEKK(",order[1],",",order[2],") is not implemented.")
	}
	
	# init the initial parameters
	params.length = count.triangular(series.count) + (order[1] * series.count^2) + (order[2] * series.count^2) 
	if(is.null(params))
	{
		if(order[1] == 1 && order[2] == 1)
		{
			if(series.count == 2)
			{
				params = c(1, 0.8, 1, 0.5, -0.4, 0, 0.3, 0.4, -0.3, 0.5, 0.8)
			}
			else
			{
				params = c(1, 0.2, 1.04, 0.3, 0.01, 0.9, 
					0.3, -0.02, -0.01, 0.01, 0.4, -0.06, 0.02, 0.3, 0.5,
					0.2, 0.01, -0.1, -0.03, 0.3, -0.06, 0.7, 0.01, 0.5)
			}
		}
		else if(order[1] == 1 && order[2] == 2)
		{
			params = c(1, 0.8, 1, 0.6, -0.4, 0, 0.4, 0.2, 0.1, -0.1, 0.3, 0.6, -0.2, 0.3, 0.5)
		}
		else if(order[1] == 2 && order[2] == 1)
		{
			params = c(1, 0.8, 1, 0.6, -0.4, 0, 0.4, 0.6, -0.2, 0.3, 0.5, 0.2, 0.1, -0.1, 0.3)
		}
		else if(order[1] == 2 && order[2] == 2)
		{
			params = c(1, 0.8, 1, 0.6, -0.4, 0, 0.3, 0.2, 0.1, -0.1, 0.3, 0.6, -0.2, 0.3, 0.5, 0.2, 0.1, -0.1, 0.3)
		}
		else if(order[1] == 0 && order[2] == 1)
		{
			params = c(1, 0.8, 1, 0.9, 0.3, -0.4, 0.8)
		}
		# fill the rest: might be a bad idea
		params = c(params, rep(0.01, params.length - length(params)))
	}
	
	# check the parameter list
	if(length(params) != params.length)
	{
		stop("length of parameters doesn't match the requiered length for the requested BEKK model");
	}
	
	# how many parameter matrices in total
	total.par.matrices = 1 + order[1] + order[2]

	# declare the parameter list
	buff.par = list()
	
	# first initialize the C matrix
	tmp.array = array(rep(0, series.count^2), dim = c(series.count, series.count))
	iter = 1
	for(i in 1:series.count)
	{
		for(j in 1:series.count)
		{
			if(i >= j)
			{
				tmp.array[j,i] = params[iter]
				iter = iter + 1
			}
		}
	}
	buff.par[[1]] = tmp.array

	# following loop initalizes the ARCH and GARCH parameter matrices respectively
	for(count in 1:(order[2] + order[1]))
	{
		buff.par[[count + 1]] = array(params[(count.triangular(series.count) + 1 + (count - 1) * series.count^2):(count.triangular(series.count) + 1 + series.count^2 + (count - 1) * series.count^2)], dim = c(series.count, series.count));
	}

	# calculate the transposes of the parameter matrices
	buff.par.transposed = list()
	for(count in 1:length(buff.par))
	{
		buff.par.transposed[[count]] = t(buff.par[[count]])
	}

	# compute the generalised kronecker product sums
	kronecker.sum = 0
	for(count in 2:total.par.matrices)
	{
		kronecker.sum = kronecker.sum + kronecker(buff.par[[count]], buff.par[[count]])
	}
	# compute the eigenvalues
	tmp.svd = svd(kronecker.sum)
	eigenvalues = tmp.svd$d

	# compute the unconditional covariance matrix
	numerat = t(buff.par[[1]]) %*% buff.par[[1]]
	dim(numerat) = c(series.count^2,1)
	denom = solve(diag(rep(1, series.count^2)) - kronecker.sum)
	sigma = denom %*% numerat
	dim(sigma) = c(series.count, series.count)

	#  simulate a two-dimensional normal white noise process:

	T = T + 50	# the first 50 are later to be discarded
	nu = rnorm(series.count * T)
	nu = array(nu, dim = c(series.count, T))

	# construct the simulated BEKK process
	HLAGS = list()
	for(count in 1:max(order))
	{
		HLAGS[[count]] = array(rep(0, series.count^2), dim = c(series.count,series.count))
		diag(HLAGS[[count]]) = 1
	}
	
	H = array(rep(0, series.count^2), dim = c(series.count, series.count))
	cor = numeric()
	eps.list = list()			# declare the estimated standard deviation series
	for(i in 1:series.count)
	{
		eps.list[[i]] = numeric()
	}
	
	sd = list()			# declare the estimated standard deviation series
	for(i in 1:series.count)
	{
		sd[[i]] = numeric()
	}

	# initialize the first instances of the time series where
	# HLAGS are not available
	for(count in 1:max(order))
	{
		for(i in 1:series.count)
		{
			eps.list[[i]][count] = 0 
		}
	}

	eps = array(rep(0, series.count), dim = c(series.count,1))
	CTERM = buff.par.transposed[[1]] %*% buff.par[[1]] # the C'C term
	for(count in (max(order) + 1):T)
	{
		# do the swap calculation for H terms
		if(order[1] >= 2)
		{
			for(tmp.count in max(order):2)
			{
				HLAGS[[tmp.count]] = HLAGS[[(tmp.count - 1)]]
			}
		}
		HLAGS[[1]] = H
		H = CTERM
		ord1 = 1
		for(tmp.count in 1:(order[2] + order[1]))
		{
			if(tmp.count <= order[2])
			{
				# ARCH EFFECT
				tmp.arr = numeric()
				for(scount in 1:series.count)
				{
					tmp.arr[scount] = eps.list[[scount]][count - tmp.count]
				}
				eps = array(tmp.arr, dim = c(series.count,1))
				H = H + buff.par.transposed[[(tmp.count + 1)]] %*% eps %*% t(eps) %*% buff.par[[(tmp.count + 1)]]
			}
			else
			{
				# GARCH EFFECT
				H = H + buff.par.transposed[[(tmp.count + 1)]] %*% HLAGS[[ord1]]  %*% buff.par[[(tmp.count + 1)]]
				ord1 = ord1 + 1
			}
		}

		svdH = svd(H)
		sqrtH = svdH$u %*% diag(sqrt(svdH$d)) %*% t(svdH$v)
		eps = sqrtH %*% nu[,count]
		#cor[count] = H[1,2]/(sqrt(H[1,1] * H[2,2]))
		for(i in 1:series.count)
		{
			sd[[i]][count] = sqrt(H[i, i]) 
			eps.list[[i]][count] = eps[i,1]
		}
	}
	
	
	names(order) <- c("GARCH component", "ARCH component")
	names(buff.par) <- as.integer(seq(1, order[1] + order[2] + 1))
	
	nu = nu[51:T]
	#cor = cor[51:T]	
	for(i in 1:series.count)
	{
		sd[[i]] = sd[[i]][51:T] 
		eps.list[[i]] = eps.list[[i]][51:T] 
	}
	T = T - 50
	
	mvBEKK.sim <- list(
		length = T,
		series.count = series.count,
		order = order,
		params = params,
		true.params = buff.par,
		eigenvalues = eigenvalues,
		uncond.cov.matrix = sigma,
		white.noise = nu,
		eps = eps.list,
		cor = cor,
		sd = sd
	)

	class(mvBEKK.sim) <- "mvBEKK.sim"
	
	cat("Class attributes are accessible through following names:\n")
	cat(names(mvBEKK.sim), "\n")
	
	return(mvBEKK.sim)
}
