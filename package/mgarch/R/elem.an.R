#####################################################################################
###   THIS PROGRAM DEFINES THE FUNCTION elem.an ("ELEMENTARY ANALYSIS").
###   IT COMPUTES DAILY AND WEEKLY RETURNS OF A TIME SERIES OF DAILY VALUES,
###   ANALYZES THE RETURN SERIES AND MAKES A LIST OF MISSING DAYS AND WEEKS.
###
###   USAGE:
###
###       elem.an(index.name, from.to = NULL, return.formula = 'simple',
###               make.bootstrap.se = TRUE, make.bull.indicator = FALSE,
###               make.weekly = TRUE, save.data.files = TRUE,
###               save.statistics = FALSE, verbose = TRUE)
###
###   ARGUMENTS:
###
###     - index.name
###     - from.to (optional): a pair of calendar dates marking the beginning and end
###       of the period to be analyzed
###     - return.formula: The method to be used to compute returns.
###       With 'simple' ('log'), simple returns in percent (log returns, respectively)
###       are computed.
###     - make.bootstrap.se: Logical. Should bootstrap standard errors of skewness
###       and kurtosis be estimated?
###     - make.bull.indicator: Logical. If TRUE, an indicator 'bull' is computed
###       which shows if a day belongs to a bull period (bull is TRUE) or to
###       a bear period (bull is FALSE)
###     - make.weekly: Should weekly returns and their statistics be computed?
###     - save.data.files: Should a file with (daily and weekly) date, weekday,
###       index value, and returns and another file with missing day/week
###       information be saved?
###     - save.statistics: Logical. If TRUE, the results of the analysis are saved
###       in a file.
###     - verbose: Logical. If TRUE, the program comments what it is doing.
###
###   DETAILS:
###
###       The function computes all returns and statistics on the basis of a
###       data frame stored in a file with name paste(index.name, '.dat', sep = '');
###       for example: dji.dat. This data frame must have the two columns 'date', with
###       dates in ISO 8601 format, and 'index', the (e.g. closing) values.
###
###       If from.to is not specified, the function computes a return series and
###       statistics for the maximum period available.
###
###       Bootstrap standard errors of skewness and kurtosis are estimated on the
###       basis of 100 bootstrap replications.
###
###       Bull indicators are obtained in a two-step procedure: Firstly, the series of
###       index values is smoothed using a one-sided moving average of length 50
###       with weights linearly decreasing in reverse time (more recent values are given
###       higher weight). Secondly, the differenced smoothed series is smoothed again with
###       a two-sided moving average of length 15. A day is said to be in a bull period
###       if this smoothed series is positive on the day in question.
###
###       The data frames missing.daily and missing.weekly also contain information
###       on the number of days (weeks) missing, and the return on the next day (week)
###       for which an index value is again available. This is meant to help the user
###       judge if a very grave event may have happened during the period of missing
###       quotes.
###
###   VALUE:
###
###   A list with components:
###     - index.name: name of index analyzed.
###     - return.formula: Method used to compute returns.
###     - from.to: Vector of two dates: beginning and end of return series.
###     - dataset.daily: Data frame with date, weekday, index and return series.
###     - missing.daily: Data frame with date, weekday, number of days missing,
###       and return with respect to the next day which is not missing
###     - statistics.daily: Data frame with statistics computed.
###     - bull.statistics: Number of days in bull/bear periods etc.
###     - temp.data.daily: Data frame; a 'calendar' with _all_ days (including
###       weekends) and sorted index values and returns.
###     - dataset.weekly: Data frame with date, weekday (Tuesday to indicate weekly
###       returns), index and return series.
###     - missing.weekly: Data frame with date, weekday, number of weeks missing,
###       and return with respect to the next week which is not missing
###     - statistics.weekly: Data frame with statistics computed.
###     - make.bootstrap.se: Was it requested to compute bootstrap standard errors
###       of skewness and kurtosis?
###     - make.bull.indicator: Was it requested to compute the bull indicator?
###     - make.weekly: Were weekly return series and statistics requested?
###     - save.data.files: Was it requested to save dataset.daily and dataset.weekly?
###     - save.statistics: Was it requested to save statistics.daily and statistics.weekly?
###     - verbose: Were comments requested?
###
###   EXAMPLES OF TYPICAL USAGE:
###
###     - e1 = elem.an('xu100')  # ----  THE FILE xu100.dat MUST EXIST!!!
###     - e1 = elem.an('xu100', from.to = c('1990-01-01', '2003-12-31'))
###     - e1 = elem.an('xu100', return.formula = 'log')
###     - e1 = elem.an('xu100', from.to = c('1990-01-01', '2003-12-31'), save.statistics = T)
###     - e1 = elem.an('xu100', make.bull.indicator = T, make.weekly = F)
###
###   Version elem.an_02.R: Uses numbers for weekdays (Monday = "1" etc.)
###   to avoid language conflicts.
###
###   (Copyright Harald Schmidbauer, 2004-09-07, 2006-04-19)
#####################################################################################


elem.an <- function(
  index.name,
  from.to = NULL, # to restrict the analysis to a subperiod
  return.formula = 'simple', # one of: 'simple', 'log'
  make.bootstrap.se = T,  # if T, computes bootstrap standard errors for skewness etc.
  make.bull.indicator = F,  # if T, makes a bull indicator for every day
  make.weekly = T,  # if T, computes weekly returns and statistics
  save.data.files = T, # if T, 'data' and 'missing' files are saved -- daily and weekly
  save.statistics = F,  # if T, statistics of the return series are saved -- daily and weekly
  verbose = T) # if T, program comments what it's doing
{
time1 = Sys.time()

################################################################################
###  CHECK IF INPUT IS OK
################################################################################
if (return.formula != 'simple' && return.formula != 'log') {
  stop("return.formula must be either 'simple' or 'log'")
}
if (is.logical(make.bootstrap.se) == F) { stop("make.bootstrap.se must be logical") }
if (is.logical(make.bull.indicator) == F) { stop("make.bull.indicator must be logical") }
if (is.logical(make.weekly) == F) { stop("make.weekly must be logical") }
if (is.logical(save.data.files) == F) { stop("save.data.files must be logical") }
if (is.logical(save.statistics) == F) { stop("save.statistics must be logical") }
if (is.logical(verbose) == F) { stop("verbose must be logical") }

################################################################################
###  DEFINE INPUT AND OUTPUT FILE NAMES:
################################################################################
index1.file = paste(index.name, '.dat', sep = '')
output.data.daily.file.name = paste('data_', index.name, '_daily.dat', sep = '')
output.statistics.daily.file.name = paste('statistics_', index.name, '_daily.dat', sep = '')
output.missing.daily.file.name = paste('missing_', index.name, '_daily.dat', sep = '')
output.data.weekly.file.name = paste('data_', index.name, '_weekly.dat', sep = '')
output.statistics.weekly.file.name = paste('statistics_', index.name, '_weekly.dat', sep = '')
output.missing.weekly.file.name = paste('missing_', index.name, '_weekly.dat', sep = '')

################################################################################
###  READ IN THE DATA, COMPUTE DAILY RETURNS:
################################################################################
dataset1 = read.table(index1.file, header = T)
if (length(dataset1$date) == 0 || length(dataset1$index) == 0 ) {
  message = "dataset not suitable -- date and/or index column missing in file"
  message = paste(message, index1.file)
  stop(message)
}

compute.simple.returns = function(series) {
  ret = diff(series)/series[1:(length(series) - 1)]*100
}

compute.log.returns = function(series) {
  ret = diff(log(series))
}

if ( return.formula == 'simple' ) { ret = compute.simple.returns(dataset1$index) }
if ( return.formula == 'log'    ) { ret = compute.log.returns(dataset1$index) }

# drop the first day of the dataset, for whch no return is available:
dataset1 = dataset1[-1,]
dataset1 = data.frame(date = dataset1$date, index = dataset1$index, ret)
if (verbose == T) {cat('Input files read, returns computed.\n')}

################################################################################
# The following function makes a dataset shorter by omitting the days
# before from.to[1] and the days after from.to[2].
################################################################################
shorten.dataset = function(from.to, dataset) {
  time.differences =
         difftime(as.POSIXct(from.to[1]), dataset$date)
  positive.time.differences = ( time.differences > 0 )
  # A positive time difference means: dataset begins earlier than from.to[1].
  # These excess days are dismissed from dataset:
  dataset = dataset[positive.time.differences == F,]
  # The same for the end:
  time.differences =
         difftime(as.POSIXct(from.to[2]), dataset$date)
  positive.time.differences = ( time.differences < 0 )
  # A negative time difference means: dataset ends later than from.to[2].
  # These excess days are dismissed from dataset:
  dataset = dataset[positive.time.differences == F,]
  dataset
}

if ( is.null(from.to) == F ) {
  dataset1 = shorten.dataset(from.to, dataset1)
}


################################################################################
###  DETERMINE THE FIRST AND THE LAST DAY AND
###  DEFINE A SEQUENCE OF DAYS FROM THE FIRST TO THE LAST DAY (ALL DAYS!!!)
################################################################################
###  1. Determine the first and the last day for which
###     both data sets are available:
first.day = as.character(as.POSIXct(dataset1$date[1]))
first.year = substring(first.day, 1, 4)
first.month = substring(first.day, 6, 7)
first.day.in.month = substring(first.day, 9, 10)

last.day = as.character(as.POSIXct(dataset1$date[length(dataset1$date)]))
last.year = substring(last.day, 1, 4)
last.month = substring(last.day, 6, 7)
last.day.in.month = substring(last.day, 9, 10)

###  2. Define the sequence:
tage = seq(ISOdate(first.year, first.month, first.day.in.month),
           ISOdate(last.year,   last.month,  last.day.in.month), 'day')
tage = trunc.POSIXt(tage, units = "days")
tage.as.character = as.character(tage, '%Y-%m-%d')
weekday = format(tage, '%w')

if ( length(tage.as.character) <= 20 ) {
  stop('Time series not long enough.', call. = F)
}

################################################################################
###  NOW MAKE A DATA FRAME (temp.data.daily!!!) WHICH CONSISTS OF ALL DATES, WEEKDAYS,
###  AND THE INDEX VALUES, SORTED BY DATES.
###  THE DATA FRAME temp.data.daily IS THE BASIS FOR ALL FURTHER DATA FRAMES.
################################################################################

###  The following function produces as output a sequence of index values
###  ordered by date, as given in tage.as.character.
###  (I.e. tage.as.character is like a calender, into which the dataset values are entered.)
sort.by.date = function(tage.as.character, dataset) {

  # Now we can begin to sort the index values, according to the 'calender'
  # defined by tage.as.character:
  days = as.character(dataset$date)
  first.year = as.numeric(substr(tage.as.character[1], 1, 4))
  last.year = as.numeric(substr(tage.as.character[length(tage.as.character)], 1, 4))
  sortedindex = numeric()
  sortedret = numeric()
  for (my.year in first.year:last.year) {
    cat('... working on year', my.year, '\n')
    # limit both series - tage.as.character and dataset - to first.year:
    year.ind.1 = as.numeric(substr(tage.as.character, 1, 4)) == my.year
    tage.as.character.temp =  tage.as.character[year.ind.1 == T]
    year.ind.2 = as.numeric(substr(days, 1, 4)) == my.year
    dataset.temp = dataset[year.ind.2 == T, ]
    days.temp = as.character(dataset.temp$date)
    sortedindex.temp =  rep(NA, length(tage.as.character.temp))
    sortedret.temp =  rep(NA, length(tage.as.character.temp))
    for (i in 1:length(days.temp)) {
      tagesindex = match(days.temp[i], tage.as.character.temp)
      sortedindex.temp[tagesindex] = dataset.temp$index[i]
      sortedret.temp[tagesindex] = dataset.temp$ret[i]
      tage.as.character.temp[tagesindex] = NA
      }
    sortedindex = c(sortedindex, sortedindex.temp)
    sortedret = c(sortedret, sortedret.temp)
  }
  array(c(sortedindex, sortedret), dim = c(length(sortedindex), 2))
}

cat('Now beginning to sort...\n')
time01 = Sys.time()
index.and.ret.sorted = sort.by.date(tage.as.character, dataset1)
time02 = Sys.time()
cat('Finished. Time spent on sorting:', difftime(time02, time01), '\n')

dataset.new = data.frame(
                date = tage.as.character,
                weekday,
                index = index.and.ret.sorted[,1],
                ret = index.and.ret.sorted[,2])
# weekend = ( dataset.new$weekday == 'Saturday' | dataset.new$weekday == 'Sunday' )
weekend = ( dataset.new$weekday == '6' | dataset.new$weekday == '0' )

dataset.new = dataset.new[weekend == F, ]
dataset.new = data.frame(
                date = dataset.new$date,
                weekday = dataset.new$weekday,
                index = dataset.new$index,
                ret = dataset.new$ret)

# which days are available:
ava = !is.na(dataset.new$index)
observations = length(ava[ava == T])
NAs = length(ava[ava == F])

if ( observations <= 10 ) {
  stop('Time series not long enough.', call. = F)
}


################################################################################
###   The following function has as argument a sequence x of:
###    - either TRUE and FALSE
###    - or 0s and 1s.
###   Output is a data frame 'analyze.iterations' with the following components:
###    - the original sequence, 'x'
###    - a sequence 'new.it' which indicates if a new iteration is beginning at a
###      certain place
###    - a sequence 'l.of.it' which contains the length of the coming iteration
###      in a place where a new iteration is beginning, and NA in other places
################################################################################

analyze.iterations = function(x) {
  new.it = ( x != c( !x[1], x[1:(length(x) - 1)]))
  new.it = ( c(x, !x[length(x)]) != c( !x[1], x))
  l.help = diff(which(new.it == T))
  l.of.it = rep(NA, length(x))
  new.it = new.it[-length(new.it)]
  l.of.it[which(new.it)] = l.help
  analyze.iterations = data.frame(x, new.it, l.of.it)
  analyze.iterations
}

result = analyze.iterations(ava)
ret.after = dataset.new[which((result$new.it) == T & (result$x == T) == T),]

missing.temp = data.frame(
                 date = dataset.new$date,
                 weekday = dataset.new$weekday,
                 days.missing = result$l.of.it,
                 ret = dataset.new$ret)

# reduce this data frame to those days which are missing:
missing.after = missing.temp[ ( result$new.it == T ) & ( result$x == T ) , ]
missing.after = missing.after[-1,]
missing.daily = missing.temp[ ( result$new.it == T ) & ( result$x == F ) , ]
missing.daily = data.frame(
              date = missing.daily$date,
              weekday = missing.daily$weekday,
              days.missing = missing.daily$days.missing,
              ret.after = missing.after$ret)

dataset.daily = data.frame(date = dataset.new$date, weekday = dataset.new$weekday,
                           index = dataset.new$index, ret = dataset.new$ret)
# temp.data.daily does not contain weekends, but all NAs.
# It is the basis for computing the weekly dataset.
temp.data.daily = dataset.daily
dataset.daily = na.omit(dataset.daily)
dataset.daily = data.frame(date = dataset.daily$date, weekday = dataset.daily$weekday,
                           index = dataset.daily$index, ret = dataset.daily$ret)
rm(dataset.new)


######################################################################
###  FOR DAILY RETURNS, WE MAKE THE BULL INDICATOR FOR THE INDEX SERIES.
######################################################################
if (make.bull.indicator == T) {
  bull = function(index) {
    fac = 1
    for (i in 2:50) {fac[i] = fac[i-1] - 0.02}
    fac = fac/sum(fac)
    f.index = filter(index, filter = fac, method = 'convolution', sides = 1)
    bull = (filter(diff(f.index), filter = rep(1/15, 15), method = 'c', sides = 2) > 0)
    bull = c(NA, bull)
  }
  bull = bull(dataset.daily$index)
  dataset.daily = data.frame(date = dataset.daily$date, weekday = dataset.daily$weekday,
                             index = dataset.daily$index, ret = dataset.daily$ret,
                             bull = bull)
  if (verbose == T) {cat('Bull indicators have been computed.\n')}
  ###  Now analyze the sequence of bull indicators:
  analyze.bull = analyze.iterations(na.omit(bull))
  bull.durations = as.numeric(na.omit(analyze.bull[analyze.bull$x == T,]$l.of.it))
  bull.days = sum(bull.durations)
  avg.bull.duration = mean(bull.durations)
  number.bull.runs = length(bull.durations)
  bear.durations = as.numeric(na.omit(analyze.bull[analyze.bull$x == F,]$l.of.it))
  bear.days = sum(bear.durations)
  avg.bear.duration = mean(bear.durations)
  number.bear.runs = length(bear.durations)
  # now produce bull statistics:
  bull.statistics = data.frame(value =
        c(bull.days, bear.days,
          formatC(avg.bull.duration, digits = 3, format = 'f'),
          formatC(avg.bear.duration, digits = 3, format = 'f'),
          number.bull.runs, number.bear.runs))

  row.names(bull.statistics) =
                        c('bull.days', 'bear.days',
                          'avg.bull.duration', 'avg.bear.duration',
                          'number.bull.runs', 'number.bear.runs')
}



################################################################################
###   THE FUNCTION make.statistics COMPUTES MEAN, VARIANCE, SKEWNESS, KURTOSIS ETC.
###   OF THE RETURNS (ret COLUMN) IN dataset AND PUTS THEM TOGETHER INTO A DATA FRAME.
################################################################################
skewness = function(x, na.rm = FALSE)
{
  if (na.rm)
      x <- x[!is.na(x)]
  mean((x - mean(x))^3)/(sd(x)^3)
}

kurtosis = function(x, na.rm = FALSE)
{
  if (na.rm)
      x <- x[!is.na(x)]
  mean((x - mean(x))^4)/(var(x)^2) - 3
}

bootstrap <- function(x, nboot, theta)
{  # simplified version of bootstrap function in package 'bootstrap'
  n <- length(x)
  bootsample <- matrix(sample(x, size = n*nboot, replace = TRUE), nrow = nboot)
  thetastar <- apply(bootsample, 1, theta)
  thetastar
}

# bootstrap estimate of mean standard error:
boot.est.se.mean = function(x, nboot)
{
  mybootstrap = bootstrap(x, nboot, mean)
  sd(na.omit(mybootstrap))
}

# bootstrap estimate of skewness standard error:
boot.est.se.skewness = function(x, nboot)
{
  mybootstrap = bootstrap(x, nboot, skewness)
  sd(na.omit(mybootstrap))
}

# bootstrap estimate of kurtosis standard error:
boot.est.se.kurtosis = function(x, nboot)
{
  mybootstrap = bootstrap(x, nboot, kurtosis)
  sd(na.omit(mybootstrap))
}

make.statistics = function(dataset)
{
  cat('Computing statistics...\n')
  nboot = 100 # this is the number of bootstrap replications
  dataset.no.NAs = na.omit(dataset)
  temp.ret = dataset.no.NAs$ret
  if ( make.bootstrap.se == T ) {
    cat('...bootstrap...')
    time.b1 = Sys.time()
    mean.se = boot.est.se.mean(temp.ret, nboot)
    skewness.se = boot.est.se.skewness(temp.ret, nboot)
    kurtosis.se = boot.est.se.kurtosis(temp.ret, nboot)
    cat('finished. Time required for bootstrap:', difftime(Sys.time(), time.b1), '\n')
  }
  if ( make.bootstrap.se == F ) {
    mean.se = sd(temp.ret)/sqrt(length(temp.ret))
    skewness.se = NA
    kurtosis.se = NA
  }
  # compute everything:
  first = as.character(dataset.no.NAs$date[1]) # first day / week
  length.dataset.no.NAs = length(dataset.no.NAs$date) # number of observations (i.e., available)
  last = as.character(dataset.no.NAs$date[length.dataset.no.NAs]) # last day / week
  NAs = length(dataset$date) - length.dataset.no.NAs   # days / weeks not available
  mean = formatC(mean(temp.ret), digits = 5, format = 'f') # mean
  mean.se = formatC(mean.se, digits = 5, format = 'f')
  var = formatC(var(temp.ret), digits = 5, format = 'f')  # variance
  sd = formatC(sd(temp.ret), digits = 5, format = 'f') # standard deviation
  skewness = formatC(skewness(temp.ret), digits = 5, format = 'f')
  if ( make.bootstrap.se == T ) { skewness.se = formatC(skewness.se, digits = 5, format = 'f') }
  if ( make.bootstrap.se == F ) { skewness.se = NA }
  kurtosis = formatC(kurtosis(temp.ret), digits = 5, format = 'f')
  if ( make.bootstrap.se == T ) { kurtosis.se = formatC(kurtosis.se, digits = 5, format = 'f')}
  if ( make.bootstrap.se == F ) { kurtosis.se = NA }
  min.ret = formatC(min(temp.ret), digits = 5, format = 'f')
  lower.quartile.ret =  formatC(quantile(temp.ret, probs = 0.25), digits = 5, format = 'f')
  median.ret =  formatC(quantile(temp.ret, probs = 0.50), digits = 5, format = 'f')
  upper.quartile.ret =  formatC(quantile(temp.ret, probs = 0.75), digits = 5, format = 'f')
  max.ret = formatC(max(temp.ret), digits = 5, format = 'f')
  day.of.min = as.character(dataset$date[which.min(dataset$ret)])
  day.of.max = as.character(dataset$date[which.max(dataset$ret)])
  cat('... computing statistics finished.\n')
  # now put them together in a vector:
  data.frame(value = c(first, last, length.dataset.no.NAs, NAs, mean, mean.se,
    var, sd, skewness, skewness.se, kurtosis, kurtosis.se,
    min.ret, lower.quartile.ret, median.ret, upper.quartile.ret, max.ret,
    day.of.min, day.of.max))
}


################################################################################
###   COMPUTE MEAN, VARIANCE, SKEWNESS, KURTOSIS ETC., PRODUCE OUTPUT:
################################################################################
statistics.daily = make.statistics(temp.data.daily)
row.names(statistics.daily) = c('first.day', 'last.day',
                          'observations', 'NAs',
                          'mean.daily', 'se.mean.daily',
                          'var.daily', 'sd.daily',
                          'skewness.daily', 'se.skewness.daily',
                          'kurtosis.daily', 'se.kurtosis.daily',
                          'min.daily', 'lower.quartile.daily', 'median.daily',
                          'upper.quartile.daily', 'max.daily',
                          'day.of.min', 'day.of.max')

if ( save.data.files == T) {
  write.table(dataset.daily, file = output.data.daily.file.name, quote = F)
  write.table(missing.daily, file = output.missing.daily.file.name, quote = F)
}
if ( save.statistics == T) {
  write.table(statistics.daily, file = output.statistics.daily.file.name, quote = F)
}

######################################################################
###  NOW WE PROCEED TO BUILDING THE SERIES OF WEEKLY VALUES.
######################################################################
if ( make.weekly == T ) {
  if ( observations <= 20 ) {
    error.message.weekly = 'Time series not long enough to compute weekly returns.\n'
    error.message.weekly = paste(error.message.weekly, 'Try make.weekly = F.')
    stop(error.message.weekly, call. = F)
  }
  if (verbose == T) {cat('Beginning to build the weekly series.\n')}
  ###  Using temp.data.daily, build temp.data.weekly successively.
  ###  First, substitute Wednesdays (or Mondays) for those Tuesdays which are missing.
  ###  The following two functions do this.
  substitute.Wednesday = function(weekday, series) {
    new.series = series
    Tuesday.missing = ( is.na(series) & (weekday == '2') )
    which.Tuesday.missing = which(Tuesday.missing == T)
    # make sure a missing Tuesday to be replaced is not the LAST day in the series:
    if ( length(Tuesday.missing[Tuesday.missing == T]) > 0 ) {
      if ( which.Tuesday.missing[length(which.Tuesday.missing)] == length(series) ) {
        which.Tuesday.missing = which.Tuesday.missing[-length(which.Tuesday.missing)]
      }
    }
    # for the missing Tuesdays, substitute a Wednesday, which is the next day:
    new.series[which.Tuesday.missing] = series[which.Tuesday.missing + 1]
    new.series
  }

  substitute.Monday = function(weekday, series) {
    new.series = series
    Tuesday.missing = ( is.na(series) & (weekday == '2') )
    which.Tuesday.missing = which(Tuesday.missing == T)
    # make sure a missing Tuesday to be replaced is not the FIRST day in the series:
    if ( length(Tuesday.missing[Tuesday.missing == T]) > 0 ) {
      if ( which.Tuesday.missing[1] == 1 ) {
        which.Tuesday.missing = which.Tuesday.missing[-1]
      }
    }
    # for the missing Tuesdays, substitute a Monday, which is the previous day:
    new.series[which.Tuesday.missing] = series[which.Tuesday.missing - 1]
    new.series
  }

  ###  Apply these two functions to the series in temp.data.daily:
  substituted.index = substitute.Wednesday(temp.data.daily$weekday, temp.data.daily$index)
  substituted.index = substitute.Monday(temp.data.daily$weekday, substituted.index)

  ###  Build a temporary data frame with weekly data:
  temp.data.weekly = data.frame(date = temp.data.daily$date,
                                weekday = temp.data.daily$weekday,
                                index1 = substituted.index)

  ###  Remove all non-Tuesdays from this data frame and build it once again
  ###  (to get the counting index straight):
  temp.data.weekly = temp.data.weekly[temp.data.weekly$weekday == '2',]
  temp.data.weekly = data.frame(date = temp.data.weekly$date,
                                weekday = temp.data.weekly$weekday,
                                index = temp.data.weekly$index1)

  ## Now we have to construct an analogue to temp.data.daily, i.e.
  ## a data frame which contains all the returns and where the NAs have
  ## not been removed yet.
  ## This is done by the following function.
  ## (The condition of its application is that there are really NAs in the dataset.)
  make.returns.in.weekly = function(dataset) {
    ## First, we omit NAs at the beginning of temp.data.weekly:
    first.available.in.temp.data.weekly =  match(F, is.na(dataset$index))
    if ( first.available.in.temp.data.weekly > 1 ) {
      dataset = dataset[-(1:first.available.in.temp.data.weekly - 1), ]
    }
    ## The rows which still have an NA are:
    rows.with.NA = which(is.na(dataset$index))
    ## Next, we fill in the remaining NAs.
    new.rows.with.NA = rows.with.NA
    index = dataset$index
    while ( length(new.rows.with.NA) >=1 ) {
      index[new.rows.with.NA] = index[new.rows.with.NA - 1]
      new.rows.with.NA = which(is.na(index))
    }
    # index is now a series without NAs. Now we can compute the returns:
    if ( return.formula == 'simple' ) { ret = compute.simple.returns(index) }
    if ( return.formula == 'log'    ) { ret = compute.log.returns(index) }

    # We substitute NAs for those returns which have been created by shifting index,
    # and which are actually zero:
    ret[rows.with.NA - 1] = NA
    dataset = dataset[-1, ]
    dataset = data.frame(date = dataset$date, weekday = dataset$weekday,
                         index = dataset$index, ret = ret)
    dataset
  }

  ###############################
  ### NOW WE HAVE IT: temp.data.weekly, which is the precise analogue to temp.data.daily:
  ###############################
  temp.data.weekly = make.returns.in.weekly(temp.data.weekly)


  # which days are available:
  ava.weekly = !is.na(temp.data.weekly$index)
  observations.weekly = length(ava.weekly[ava.weekly == T])
  NAs.weekly = length(ava.weekly[ava.weekly == F])

  result.weekly = analyze.iterations(ava.weekly)

  missing.weekly.temp = data.frame(
                   date = temp.data.weekly$date,
                   weekday = temp.data.weekly$weekday,
                   weeks.missing = result.weekly$l.of.it,
                   ret = temp.data.weekly$ret)

  missing.weekly.after =
     missing.weekly.temp[ ( result.weekly$new.it == T ) & ( result.weekly$x == T ) , ]
  missing.weekly.after = missing.weekly.after[-1,]
  ret.after = missing.weekly.after$ret

  missing.weekly =
     missing.weekly.temp[ ( result.weekly$new.it == T ) & ( result.weekly$x == F ) , ]
  # the last return value in missing.weekly.after might be missing. In this case:
  if ( length(ret.after) < length(missing.weekly$ret) ) {
    ret.after = c(ret.after, NA)
  }

  missing.weekly = data.frame(
                date = missing.weekly$date,
                weekday = missing.weekly$weekday,
                weeks.missing = missing.weekly$weeks.missing,
                ret.after = ret.after)


  dataset.weekly = na.omit(temp.data.weekly)
  dataset.weekly = dataset.weekly[-1, ]

  dataset.weekly = data.frame(date = dataset.weekly$date, weekday = dataset.weekly$weekday,
                             index = dataset.weekly$index, ret = dataset.weekly$ret)

  statistics.weekly = make.statistics(temp.data.weekly)
  row.names(statistics.weekly) = c('first.week', 'last.week',
                          'observations', 'NAs',
                          'mean.weekly', 'se.mean.weekly',
                          'var.weekly', 'sd.weekly',
                          'skewness.weekly', 'se.skewness.weekly',
                          'kurtosis.weekly', 'se.kurtosis.weekly',
                          'min.weekly', 'lower.quartile.weekly', 'median.weekly',
                          'upper.quartile.weekly', 'max.weekly',
                          'week.of.min', 'week.of.max')

  if ( save.data.files == T) {
    write.table(dataset.weekly, file = output.data.weekly.file.name, quote = F)
    write.table(missing.weekly, file = output.missing.weekly.file.name, quote = F)
  }
  if ( save.statistics == T) {
    write.table(statistics.weekly, file = output.statistics.weekly.file.name, quote = F)
  }
}  # end of weekly stuff

if ( make.weekly == F ) {
  dataset.weekly = NULL
  missing.weekly = NULL
  statistics.weekly = NULL
}

if ( make.bull.indicator == F ) {
  bull.statistics = NULL
}


cat('Total time required: ', difftime(Sys.time(), time1), '\n')

elem.an <- list(
  index.name = index.name,
  return.formula = return.formula,
  from.to = from.to,
  dataset.daily = dataset.daily,
  missing.daily = missing.daily,
  statistics.daily = statistics.daily,
  bull.statistics = bull.statistics,
  temp.data.daily = temp.data.daily,
  dataset.weekly = dataset.weekly,
  missing.weekly = missing.weekly,
  statistics.weekly = statistics.weekly,
  make.bootstrap.se = make.bootstrap.se,
  make.bull.indicator = make.bull.indicator,
  make.weekly = make.weekly,
  save.data.files = save.data.files,
  save.statistics = save.statistics,
  verbose = verbose)

class(elem.an) = "elem.an"
cat("Class attributes can be called via the following names:\n")
cat(names(elem.an), "\n")
return(elem.an)

} # end of definition of function elem.an
