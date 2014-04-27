\name{mvBEKK.diag}

\alias{mvBEKK.diag}

\title{Support diagnosis of BEKK(p,q) model fitting}

\description{
	mvBEKK.diag prints the results of an estimation of a BEKK(p,q) model
	in a fancy format
}

\usage{
mvBEKK.diag(estimation)
}

\arguments{
  \item{estimation}{estimation object returned from the mvBEKK.est routine}
}

\details{
  \code{mvBEKK.diag} prints the results of an estimation of a BEKK(p,q) model
	in a fancy format
}

\value{
	NONE
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
  \code{\link{mvBEKK.est}} for estimation of BEKK models.
}

\examples{
\dontrun{
sim = mvBEKK.sim(series.count = 3, T = 1000) # simulate a 3 dimensional mgarch model with length of 1000
eps = data.frame(sim$eps[[1]], sim$eps[[2]], sim$eps[[3]]) # encapsulate
est = mvBEKK.est(eps) # estimate the simulated model
mvBEKK.diag(est) # print diagnosis
}}

\keyword{ts}
