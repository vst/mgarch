\name{mGJR.est}

\alias{mGJR.est}

\title{Estimate a mGJR(p,q,g) model}

\description{
	mGJR.est estimates a mGJR(p,q,g) model
	for two given time series
}

\usage{
mGJR.est(eps1, eps2, order = c(1,1,1), params = NULL, fixed = NULL, method = "BFGS")
}

\arguments{
  \item{eps1}{first time series}
  \item{eps2}{second time series}
  \item{order}{a three dimensional integer vector giving the orders of the
    model to be fitted. \code{order[2]} refers to the ARCH order and
    \code{order[1]} to the GARCH order and \code{order[3]} to the GJR order.}
  \item{params}{initial parameter list for the optimization process}
  \item{fixed}{a two dimensional vector that contains the user specified fixed parameter values.}
  \item{method}{method to be used in the optimization process. See ?optim for available options.}
}

\details{
  \code{mGJR.est} estimates a mGJR(p,q) model, where p stands for the 
  GARCH order, and q stands for the ARCH order
}

\value{
  A list of class \code{"mGJR.est"} with the following elements:
	\item{eps1}{first time series}
	\item{eps2}{second time series}
	\item{length}{length of each series}
	\item{order}{order of the mGJR model fitted}
	\item{estimation.time}{time to complete the estimation process}
	\item{total.time}{time to complete the whole routine within the mGJR.est process}
	\item{estimation}{estimation object returned from the optimization process, using \code{optim}}
	\item{aic}{the AIC value of the fitted model}
	\item{est.params}{estimated parameter matrices}
	\item{asy.se.coef}{asymptotic theory estimates of standard errors of estimated parameters}
	\item{cor}{estimated conditional correlation series}
	\item{sd1}{first estimated conditional standard deviation series}
	\item{sd2}{second estimated conditional standard deviation series}
	\item{H.estimated}{estimated series of covariance matrices}
	\item{eigenvalues}{estimated eigenvalues for sum of Kronecker products}
	\item{uncond.cov.matrix}{estimated unconditional covariance matrix}
	\item{resid1}{first estimated series of residuals}
	\item{resid2}{second estimated series of residuals}
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
}

\examples{
\dontrun{
sim = BEKK.sim(1000)
est = mGJR.est(sim$eps1, sim$eps2)
}}
\keyword{ts}
