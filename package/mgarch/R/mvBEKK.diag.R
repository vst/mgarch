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

mvBEKK.diag <-
function(
			estimation		# result returned from the BEKK.estimation process
		)
{
	cat("\tNumber of estimated series : ", length(estimation$eps),    "\n")	
	cat("\tLength of estimated series : ", estimation$series.length,    "\n")	
	cat("\tEstimation Time            : ", estimation$estimation.time,  "\n")	
	cat("\tTotal Time                 : ", estimation$total.time,       "\n")	
	cat("\tBEKK order                 : ", estimation$order,            "\n")	
	cat("\tEigenvalues                : ", estimation$eigenvalues,      "\n")	
	cat("\taic                        : ", estimation$aic,              "\n")
	cat("\tunconditional cov. matrix  : ", estimation$uncond.cov.mat,   "\n")	
	
	for(i in 1:length(estimation$eps))
	{
		cat("\tvar(resid", i, ")                : ", var(estimation$residuals[[i]]),      "\n")
		cat("\tmean(resid", i, ")               : ", mean(estimation$residuals[[i]]),     "\n")
	}
	#cat("\tcor(resid1, resid2)        : ", cor(estimation$resid1, estimation$resid2), "\n")
	
	cat("\tEstimated parameters       :\n\n")
	cat("\tC estimates:\n")
	print(estimation$est.params[[1]])
	
	if(estimation$order[2] > 0)
	{
		cat("\n\tARCH estimates:\n")
		for(count in 1:estimation$order[2])
		{
			print(estimation$est.params[[count + 1]])
		}
	}
	else
	{
		count = 0
	}
	
	if(estimation$order[1] > 0)
	{
		cat("\n\tGARCH estimates:\n")
		for(count2 in 1:estimation$order[1])
		{
			print(estimation$est.params[[(count + 1) + count2]])
		}
	}

	cat("\n\tasy.se.coef                : \n\n")
	cat("\tC estimates, standard errors:\n")
	print(estimation$asy.se.coef[[1]])
	
	if(estimation$order[2] > 0)
	{
		cat("\n\tARCH estimates, standard errors:\n")
		for(count in 1:estimation$order[2])
		{
			print(estimation$asy.se.coef[[count + 1]])
		}
	}
	else
	{
		count = 0
	}
	
	if(estimation$order[1] > 0)
	{
		cat("\n\tGARCH estimates, standard errors:\n")
		for(count2 in 1:estimation$order[1])
		{
			print(estimation$asy.se.coef[[(count + 1) + count2]])
		}
	}
	
	browser()
	
#	plot(
#			min(min(estimation$resid1),min(estimation$resid2)):max(max(estimation$resid1),max(estimation$resid2)),
#			min(min(estimation$resid1),min(estimation$resid2)):max(max(estimation$resid1),max(estimation$resid2)),
#			type = "n",
#			xlab = "resid1",
#			ylab = "resid2"
#		)
#	
#	points(estimation$resid1, estimation$resid2, pch = 21)
	
	for(i in 1:length(estimation$eps))
	{
		plot(estimation$residuals[[i]])
		browser()
		dev.off()
	}
}
