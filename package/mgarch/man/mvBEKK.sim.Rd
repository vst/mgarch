\name{mvBEKK.sim}

\alias{mvBEKK.sim}

\title{Simulate a BEKK(p,q) model}

\description{
	mvBEKK.sim simulates an N dimensional BEKK(p,q) model
	for the given length, order list, and initial
	parameter list where N is also specified by the user.
}

\usage{
mvBEKK.sim(series.count, T, order = c(1,1), params = NULL)
}

\arguments{
  \item{series.count}{number of the series to be simulated.}
  \item{T}{length of the series to be simulated.}
  \item{order}{an integer vector of length 2 giving the orders of the
    model to fit.  \code{order[2]} refers to the ARCH order and
    \code{order[1]} to the GARCH order.}
  \item{params}{a vector containing a sequence of parameter matrices' values}
}

\details{
  \code{mvBEKK.sim} simulates a BEKK(p,q) model, where p stands for the 
  GARCH order and q stands for the ARCH order
}

\value{
  A list of class \code{"mvBEKK.sim"} with the following elements:
	\item{length}{length of the series simulated}
	\item{order}{order of the BEKK model}
	\item{params}{a vector of the selected parameters}
	\item{true.params}{list of parameters in matrix form}
	\item{eigenvalues}{computed eigenvalues for sum of Kronecker products}
	\item{uncond.cov.matrix}{unconditional covariance matrix of the process}
	\item{white.noise}{white noise series used for simulating the process}
	\item{eps}{a list of simulated series}
	\item{cor}{list of series of conditional correlations}
	\item{sd}{list of series of conditional standard deviations}
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
  \code{\link{mvBEKK.est}} for estimation of BEKK models;
  \code{\link{mvBEKK.diag}} for printing out mvBEKK.est in a fancy format.
}

\examples{
\dontrun{
sim = mvBEKK.sim(series.count = 3, T = 2500)
}}

\keyword{ts}
