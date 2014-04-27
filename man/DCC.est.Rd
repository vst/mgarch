\name{DCC.est}
\alias{DCC.est}
\title{Estimates a DCC-GARCH model.}
\usage{
DCC.est(eps, params = NULL, M = 5, verbose = F)
}
\arguments{
  \item{eps}{time series}

  \item{params}{initial parameters for the optim function}

  \item{M}{the M parameter as it appears in the paper}

  \item{verbose}{predicate which defines if we go verbose
  debugging}
}
\value{
an instance of \code{DCC.est} classx
}
\description{
Estimates a DCC-GARCH model.
}

