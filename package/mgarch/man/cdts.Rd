\name{cdts}

\alias{cdts}

\title{Combine daily time series}

\description{
	cdts combines daily time series in a daily or weekly returns format. 
}

\usage{
cdts(file.names, return.formula = 'simple', from.to, daily.availability = 1, weekly.availability = 1, verbose = T)
}

\arguments{
  \item{file.names}{a vector containing the level series (``.dat'' extension should be omitted)}
  \item{return.formula}{return formula used: 'simple' or ???}
  \item{from.to}{a vector for specifying the date range. Defaults to the
  first and last available dates}
  \item{daily.availability}{indicates whether daily combination should take place or not}
  \item{weekly.availability}{indicates whether weekly combination should take place or not}
  \item{verbose}{boolean value that tells whether verbose output should be printed or not} 
}

\details{
  \code{cdts} combines daily time series in a daily or weekly returns format. It takes the name of the files (``.dat'' extension should be omitted) of the level series and combines them by calculating the returns.
}

\value{
  \item{file.names}{the name of the input files (extension omitted)} 
  \item{data.daily}{daily data} 
  \item{data.weekly}{weekly data} 
  \item{all.data.daily}{all daily data} 
}

\author{
Harald SCHMIDBAUER \email{harald@hs-stat.com}
}

\examples{
	\dontrun{combined = cdts(c("xu100", "dji", "gdaxi"))}
}

\keyword{ts}
