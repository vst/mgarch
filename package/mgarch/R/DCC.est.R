## Copyright (C) 2009 Harald SCHMIDBAUER, Angi Roesch, Vehbi Sinan
## TUNALIOGLU
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

DCC.est <-
function(
         eps,              # time series as a 2 dimensional data frame
         params = NULL,    # initial parameters for the optim function
         M = 5,            # M as it appears in the Paper
         verbose = F
	)
{
  
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

  if(series.count != 2){
    stop("There should be exactly 2 timeseries")
  }
  

  if(is.null(params)) {
    # TODO
    # these are meaningless parameters.
    # set some useful initial parameters.
    params = c(0.2,0.7)
    out("\nWarning: initial values for the parameters are set to:\n\t", params,"\n")
  }

  if(length(params) != 2){
    stop("There should be exactly 2 parameters")
  }
  
	

  # begin estimation process
  # first log the start time
  start = Sys.time()		
  out("* Starting estimation process.\n")

  est.DCCmgarch = function(M,ini.params){
    # objective function
	enc.DCCmgarchloglike <- function(theta) {
      theta1 = theta[1]
      theta2 = theta[2]
      
      Psi = array(c(1,0,0,1), dim = c(2,2,T))
      H   = array(c(0,0,0,0), dim = c(2,2,T))
      R   = array(c(1,0,0,1), dim = c(2,2,T))
      
      L = 0
      for ( i in 2:(M+1) ) { # here, we have no Psi yet! 
        L = L + as.numeric(dmvnorm(eps[i,], mean = c(0,0), sigma = diag(c(h1[i], h2[i])), log = T))
      }
      
      for ( i in (M+2):T ) { # Now, we can also use Psi!
        Psi[1,2,i-1] = cor(u1[(i-M):(i-1)], u2[(i-M):(i-1)])
        Psi[2,1,i-1] = Psi[1,2,i-1]
        R[,,i] = (1 - theta1 - theta2) * cor(eps) + theta1 * Psi[,,i-1] + theta2 * R[,,i-1]
        D = diag(sqrt(c(h1[i], h2[i])))
        H[,,i] = D %*% R[,,i] %*% D
        L = L + as.numeric(dmvnorm(eps[i,], mean = c(0,0), sigma = H[,,i], log = T))
      }
      return(-L)
    }

    # gradient
    grad.DCCmgarchloglike <- function(theta) {
      theta1 = theta[1]
      theta2 = theta[2]
      
      # initialization
      grad1 = 0
      grad2 = 0
      Psi = array(c(1,0,0,1), dim = c(2,2,T))
      H   = array(c(0,0,0,0), dim = c(2,2,T))
      R   = array(c(1,0,0,1), dim = c(2,2,T))
      dR1 = array(c(0,0,0,0), dim = c(2,2))
      dR2 = array(c(0,0,0,0), dim = c(2,2))
      
      for ( i in (M+2):T ) { # Here, we can calculate Psi!
		Psi[1,2,i-1] = cor(u1[(i-M):(i-1)], u2[(i-M):(i-1)])
		Psi[2,1,i-1] = Psi[1,2,i-1]
		R[,,i] = (1 - theta1 - theta2) * cor(eps) + theta1 * Psi[,,i-1] + theta2 * R[,,i-1]
      }
      
      for ( i in 2:T ) {
        u = array(c(u1[i],u2[i]), dim = c(2,1))
        R.inverse = solve(R[,,i])
		dR1  = -cor(eps) + Psi[,,i-1] + theta2 * dR1
		dR2  = -cor(eps) + R[,,i-1] + theta2 * dR2
		grad1 = grad1 + sum(diag(R.inverse %*% dR1)) 
		grad1 = grad1 - t(u) %*% R.inverse %*% dR1 %*% R.inverse %*% u
		grad2 = grad2 + sum(diag(R.inverse %*% dR2))
		grad2 = grad2 - t(u) %*% R.inverse %*% dR2 %*% R.inverse %*% u
      }
      return(c(grad1/2,grad2/2))
	}
    
	# constraint optimization of DCC parameters
	result =
      constrOptim(ini.params,
                  enc.DCCmgarchloglike,
                  grad.DCCmgarchloglike,
                  ui= rbind(c(1,0),c(0,1),c(-1,-1)),
                  ci = c(0,0,-1))
	
    # one iteration of optim, just to get the hessian returned (which is not provided by
	# constrOptim):
	resulth =
      optim(result$par,
            enc.DCCmgarchloglike,
            grad.DCCmgarchloglike,
            control = list(maxit=1),
            hessian = TRUE)

	if (is.finite(det(resulth$hessian))) {
      Sigma = solve(resulth$hessian)
      se = c(sqrt(diag(Sigma)))
	}
	else {
      se = rep('NaN',2)
	}
    est.DCCmgarch = list(par = result$par, se = se)
  }

  T = length(eps[,1])
  ini.params = params

  out('\n#########################################################\n')
  out('DCC MGarch estimation, model by Tse and Tsui (2002)')
  out('\n#########################################################\n\n')
  out('Number of time series:',length(eps[1,]),'\n')
  out('Length of time series:',T,'\n')
  out('Model specification: M=',M,'\n')
  out('Initial DCC parameter values:',ini.params,'\n\n')


  ############################################################################ 
  # First step in the estimation: fit a GARCH to the individual time series
  ############################################################################ 
  
  my.garch1 = garch(eps[,1], start = c(1,0.5,0.5), trace = F)
  my.garch2 = garch(eps[,2], start = c(1,0.5,0.5), trace = F)
  # define estimated parameters:
  a10 = as.numeric(my.garch1$coef[1])
  a11 = as.numeric(my.garch1$coef[2])
  b11  = as.numeric(my.garch1$coef[3])
  a20 = as.numeric(my.garch2$coef[1])
  a21 = as.numeric(my.garch2$coef[2])
  b21  = as.numeric(my.garch2$coef[3])
  
  out('\n#########################################################\n')
  out('Univariate Garch parameter estimates:')
  out('\n#########################################################\n\n')
  out('Series 1:\n')
  out(summary(my.garch1))
  
  out('#########################################################\n\n')

  out('Series 2:\n')
  out(summary(my.garch2))

  # series of fitted conditional variances:
  h1 = my.garch1$fitted.values[,1]^2
  h2 = my.garch2$fitted.values[,1]^2
  
  u1 = eps[,1] / sqrt(h1)
  u2 = eps[,2] / sqrt(h2)

  ##############################################################################
  # Second step in the estimation: estimation of theta1 and theta2
  ##############################################################################
  
  out('#########################################################\n\n')
  out('DCC Garch parameter estimates:\n')
  out('Now working on DCC parameter estimation...\n')

  result.est.DCCmgarch = est.DCCmgarch(M,ini.params)
  out(result.est.DCCmgarch)

  out('#########################################################\n\n')

  # estimation completed
  out("* Estimation process completed.\n")

  # log estimation time
  est.time = difftime(Sys.time(), start)


  DCC.est <- list(
                     eps = eps,
                     series.length = series.length,
                     estimation.time = est.time,
                     total.time = difftime(Sys.time(), start),
                     result = result.est.DCCmgarch
                     )

  class(DCC.est) = "DCC.est"
  
  out("Class attributes are accessible through following names:\n")
  out(names(DCC.est), "\n")
  
  return(DCC.est)
}
