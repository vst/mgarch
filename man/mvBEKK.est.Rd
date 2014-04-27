\name{mvBEKK.est}

\alias{mvBEKK.est}

\title{Estimate a BEKK(p,q) model}

\description{
	mvBEKK.est estimates a BEKK(p,q) model
	for two given time series
}

\usage{
mvBEKK.est(eps, order = c(1,1), params = NULL, fixed = NULL, method = "BFGS", verbose = F)
}

\arguments{
  \item{eps}{a data frame contaning all the time series}
  \item{order}{an integer vector of length 2 giving the orders of the
    model to be fitted. \code{order[2]} refers to the ARCH order and
    \code{order[1]} to the GARCH order.}
  \item{params}{a vector of initial parameters for the optimization process}
  \item{fixed}{a two dimensional array that contains the user specified fixed parameter values.}
  \item{method}{method to be used in the optimization process. See \code{optim} function for available options.}
  \item{verbose}{If set to True, verbose output will be printed on the screen}
}

\details{
  \code{mvBEKK.est} estimates a BEKK(p,q) model, where p stands for the 
  GARCH order, and q stands for the ARCH order
}

\value{
  A list of class \code{"mvBEKK.est"} with the following elements:
  	\item{eps}{a data frame contaning all time series}
	\item{length}{length of the series}
	\item{order}{order of the BEKK model fitted}
	\item{estimation.time}{time to complete the estimation process}
	\item{total.time}{time to complete the whole routine within the mvBEKK.est process}
	\item{estimation}{estimation object returned from the optimization process, using \code{optim}}
	\item{aic}{the AIC value of the fitted model}
	\item{est.params}{list of estimated parameter matrices}
	\item{asy.se.coef}{list of asymptotic theory estimates of standard errors of estimated parameters}
	\item{cor}{list of estimated conditional correlation series}
	\item{sd}{list of estimated conditional standard deviation series}
	\item{H.estimated}{list of estimated series of covariance matrices}
	\item{eigenvalues}{estimated eigenvalues for sum of Kronecker products}
	\item{uncond.cov.matrix}{estimated unconditional covariance matrix}
	\item{residuals}{list of estimated series of residuals}
}

\references{
	Bauwens L., S. Laurent, J.V.K. Rombouts:
	Multivariate GARCH models: A survey, April, 2003
	
	Bollerslev T.:
	Modelling the coherence in short-run nominal exchange rate: A multivariate generalized ARCH approach, Review of Economics and Statistics, 498--505, 72, 1990
	
	Engle R.F., K.F. Kroner:
	Multivariate simultaneous generalized ARCH, Econometric Theory, 122-150, 1995
	
	Engle R.F.:
	Dynamic conditional correlation: A new simple class of multivariate GARCH models, Journal of Business and Economic Statistics, 339--350, 20, 2002

	Tse Y.K., A.K.C. Tsui:
	A multivariate generalized autoregressive conditional heteroscedasticity model with time-varying correlations, Journal of Business and Economic Statistics, 351-362, 20, 2002
}

\author{
Harald SCHMIDBAUER \email{harald@hs-stat.com},
Vehbi Sinan TUNALIOGLU \email{vst@vsthost.com}
}

\seealso{
  \code{\link{mvBEKK.sim}} for simulation of BEKK models;
  \code{\link{mvBEKK.diag}} for printing out mvBEKK.est estimation object in a fancy format.
}

\examples{
\dontrun{
sim = mvBEKK.sim(series.count = 3, T = 1000) # simulate a 3 dimensional mgarch model with length of 1000
eps = data.frame(sim$eps[[1]], sim$eps[[2]], sim$eps[[3]]) # encapsulate
est = mvBEKK.est(eps) # estimate the simulated model
}}

\keyword{ts}
