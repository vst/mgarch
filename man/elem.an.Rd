\name{elem.an}

\alias{elem.an}

\title{Elementary analysis}

\description{
This program defines the function elem.an ("elementary analysis").  It computes daily and weekly returns of a time series of daily values, analyzes the return series and makes a list of missing days and weeks.
}

\usage{
	elem.an(index.name, from.to = NULL, return.formula = 'simple',
		make.bootstrap.se = T, make.bull.indicator = F,
		make.weekly = T, save.data.files = T,
		save.statistics = F, verbose = T)
}

\arguments{
	\item{index.name}{name of the .dat file.}
	\item{from.to}{(optional) a pair of calendar dates marking the beginning and end of the period to be analyzed}
	\item{return.formula}{The method to be used to compute returns.  With 'simple' ('log'), simple returns in percent (log returns, respectively) are computed.}
	\item{make.bootstrap.se}{Logical. Should bootstrap standard errors of skewness and kurtosis be estimated?}
	\item{make.bull.indicator}{Logical. If TRUE, an indicator 'bull' is computed which shows if a day belongs to a bull period (bull is TRUE) or to a bear period (bull is FALSE)}
	\item{make.weekly}{Should weekly returns and their statistics be computed?}
	\item{save.data.files}{Should a file with (daily and weekly) date, weekday, index value, and returns and another file with missing day/week information be saved?}
	\item{save.statistics}{Logical. If TRUE, the results of the analysis are saved in a file.}
	\item{verbose}{Logical. If TRUE, the program comments what it is doing.}
}

\details{
  \code{elem.an}
	The function computes all returns and statistics on the basis of a
	data frame stored in a file with name paste(index.name, '.dat', sep = '');
	for example: dji.dat. This data frame must have the two columns 'date', with
	dates in ISO 8601 format, and 'index', the (e.g. closing) values.

	If from.to is not specified, the function computes a return series and
	statistics for the maximum period available.

	Bootstrap standard errors of skewness and kurtosis are estimated on the
	basis of 100 bootstrap replications.

	Bull indicators are obtained in a two-step procedure: Firstly, the series of
	index values is smoothed using a one-sided moving average of length 50
	with weights linearly decreasing in reverse time (more recent values are given
	higher weight). Secondly, the differenced smoothed series is smoothed again with
	a two-sided moving average of length 15. A day is said to be in a bull period
	if this smoothed series is positive on the day in question.

	The data frames missing.daily and missing.weekly also contain information
	on the number of days (weeks) missing, and the return on the next day (week)
	for which an index value is again available. This is meant to help the user
	judge if a very grave event may have happened during the period of missing
	quotes.
  
}

\value{
 \item{index.name}{name of index analyzed.}
 \item{return.formula}{Method used to compute returns.}
 \item{from.to}{Vector of two dates: beginning and end of return series.}
 \item{dataset.daily}{Data frame with date, weekday, index and return series.}
 \item{missing.daily}{Data frame with date, weekday, number of days missing, and return with respect to the next day which is not missing}
 \item{statistics.daily}{Data frame with statistics computed.}
 \item{bull.statistics}{Number of days in bull/bear periods etc.}
 \item{temp.data.daily}{Data frame; a 'calendar' with \_all\_ days (including weekends) and sorted index values and returns.}
 \item{dataset.weekly}{Data frame with date, weekday (Tuesday to indicate weekly returns), index and return series.}
 \item{missing.weekly}{Data frame with date, weekday, number of weeks missing, and return with respect to the next week which is not missing}
 \item{statistics.weekly}{Data frame with statistics computed.}
 \item{make.bootstrap.se}{Was it requested to compute bootstrap standard errors of skewness and kurtosis?}
 \item{make.bull.indicator}{Was it requested to compute the bull indicator?}
 \item{make.weekly}{Were weekly return series and statistics requested?}
 \item{save.data.files}{Was it requested to save dataset.daily and dataset.weekly?}
 \item{save.statistics}{Was it requested to save statistics.daily and statistics.weekly?}
 \item{verbose}{Were comments requested?}
}

\author{
Harald SCHMIDBAUER \email{harald@hs-stat.com}
}

\examples{
	\dontrun{
		e1 = elem.an('xu100')  # ----  THE FILE xu100.dat MUST EXIST!!!
		e2 = elem.an('xu100', from.to = c('1990-01-01', '2003-12-31'))
		e3 = elem.an('xu100', return.formula = 'log')
		e4 = elem.an('xu100', from.to = c('1990-01-01', '2003-12-31'), save.statistics = T)
		e5 = elem.an('xu100', make.bull.indicator = T, make.weekly = F)
	}
}

\keyword{ts}
