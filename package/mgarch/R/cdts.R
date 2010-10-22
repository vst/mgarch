cdts <- function(
               file.names,
               return.formula = 'simple',
               from.to,
               daily.availability = 1, # 1 = at least one must be available
               weekly.availability = 1, # 1 = at least one must be available
               verbose = T )
{

start.time = Sys.time()

n.files = length(file.names)

# read in the data, construct a list:
all.data = list()

for ( i in 1:length(file.names)) {
  all.data[[file.names[i]]] =
      read.table(paste(file.names[i], '.dat', sep = ''), header = T)
}

# browser()

###  1. Determine the first and the last day for which
###     all data sets are available:
begin.day = as.POSIXct(all.data[[1]]$date[1])
for ( i in 1:n.files ) {
  if (difftime(all.data[[i]]$date[1], begin.day) > 0 ) {
    begin.day = all.data[[i]]$date[1]
  }
}
# begin.day is the later of the two first days first.day.dataset1, first.day.dataset2!
begin.day = as.character(begin.day)
begin.year = substring(begin.day, 1, 4)
begin.month = substring(begin.day, 6, 7)
begin.day.in.month = substring(begin.day, 9, 10)
last.day = as.POSIXct(all.data[[1]]$date[length(all.data[[1]][,2])])
for ( i in 1:n.files ) {
  if (difftime(all.data[[i]]$date[length(all.data[[i]][,2])], last.day) < 0 ) {
    last.day = all.data[[i]]$date[length(all.data[[i]][,2])]
  }
}
# last.day is the earlier of the two last days last.day.dataset1, last.day.dataset2!
last.day = as.character(last.day)
last.year = substring(last.day, 1, 4)
last.month = substring(last.day, 6, 7)
last.day.in.month = substring(last.day, 9, 10)


###  2. Define the sequence:
tage = seq(ISOdate(begin.year, begin.month, begin.day.in.month),
           ISOdate(last.year,   last.month,  last.day.in.month), 'day')
tage = trunc.POSIXt(tage, units = "days")
tage.as.character = as.character(tage, '%Y-%m-%d')
weekday = format(tage, '%w')  # weekday = weekdays(tage)



# browser()

############################################################################################
############################################################################################
############################################################################################

# from.to = c('2000-02-05', '2000-04-26')

shorten.dataset = function(from.to, dataset) {
  time.differences =
         difftime(as.POSIXct(from.to[1]), as.POSIXct(dataset$date))
  positive.time.differences = ( time.differences > 0 )
  # A positive time difference means: dataset begins earlier than from.to[1].
  # These excess days are dismissed from dataset:
  dataset = dataset[positive.time.differences == F,]
  # The same for the end:
  time.differences =
         difftime(as.POSIXct(from.to[2]), as.POSIXct(dataset$date))
  positive.time.differences = ( time.differences < 0 )
  # A negative time difference means: dataset ends later than from.to[2].
  # These excess days are dismissed from dataset:
  dataset = dataset[positive.time.differences == F,]
  dataset
}

from.to = c(tage.as.character[1], tage.as.character[length(tage.as.character)])

all.data.short = list()
for ( i in 1:length(file.names)) {
  all.data.short[[file.names[i]]] = shorten.dataset(from.to, all.data[[i]])
}

all.data = all.data.short

# browser()

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
  for (my.year in first.year:last.year) {
    cat('... working on year', my.year, '\n')
    # limit both series - tage.as.character and dataset - to first.year:
    year.ind.1 = as.numeric(substr(tage.as.character, 1, 4)) == my.year
    tage.as.character.temp =  tage.as.character[year.ind.1 == T]
    year.ind.2 = as.numeric(substr(days, 1, 4)) == my.year
    dataset.temp = dataset[year.ind.2 == T, ]
    days.temp = as.character(dataset.temp$date)
    sortedindex.temp =  rep(NA, length(tage.as.character.temp))
    for (i in 1:length(days.temp)) {
      tagesindex = match(days.temp[i], tage.as.character.temp)
      sortedindex.temp[tagesindex] = dataset.temp$index[i]
      tage.as.character.temp[tagesindex] = NA
      }
    sortedindex = c(sortedindex, sortedindex.temp)
  }
  sortedindex
}

temp.data.daily.list = list()
temp.data.daily.list[['date']] = tage.as.character
temp.data.daily.list[['weekday']] = weekday
for ( i in 1:length(file.names)) {
  cat()
  temp.data.daily.list[[file.names[i]]] = sort.by.date(tage.as.character, all.data[[i]])
}

cat('DEBUG: temp.data.daily.list ready\n')
# browser()


################################################################################
###  NEXT, ON THE BASIS OF THE LIST temp.data.daily.list,
###  MAKE A LIST WHICH CONTAINS THE DAILY VALUES.
################################################################################

###  depending on daily.availability, decide which days have to be dismissed:
# day.available.1 is TRUE if and only if at least one day's index is available:
day.available.1 = rep(FALSE, length(temp.data.daily.list[[3]]))
for ( i in 1:length(file.names) ) {
  day.available.1 = (day.available.1 | !is.na(temp.data.daily.list[[2 + i]]))
}
# day.available.2 is TRUE if and only if all day's indices are available:
day.available.2 = rep(TRUE, length(temp.data.daily.list[[3]]))
for ( i in 1:length(file.names) ) {
  day.available.2 = (day.available.2 & !is.na(temp.data.daily.list[[2 + i]]))
}

# In any case, we dismiss all days before the first day for which BOTH indices are available:
# available.from = min(which(day.available.2 == T))
if (daily.availability == 1) {day.available = day.available.1}
if (daily.availability == 2) {day.available = day.available.2}
# day.available =
#    c(rep(FALSE, available.from - 1), day.available[available.from:length(day.available)])
if (verbose == T) {cat('Decision which days to be dismissed has been made.\n')}

###  Record which days are not available:
################################################################################
###  TO BE DONE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###
################################################################################

# browser()

###  Now dismiss those days for which no data are available:
reduced.temp.data.daily.list = list()
reduced.temp.data.daily.list[['date']] = temp.data.daily.list[['date']][day.available == T]
reduced.temp.data.daily.list[['weekday']] = temp.data.daily.list[['weekday']][day.available == T]
for ( i in 1:length(file.names)) {
  reduced.temp.data.daily.list[[file.names[i]]] =
    temp.data.daily.list[[file.names[i]]][day.available == T]
}
if (verbose == T) {cat('Not available days have been dismissed.\n')}

# browser()

###  If daily.availability == 1:
###  Shift indices to fill remaining NAs (where at least one index is available),
###  and modify the data frame reduced.temp.data.daily.
###  (This works because index.series.with.NAs cannot begin with an NA!)
eliminate.NAs = function(index.series.with.NAs) {
  index.series.without.NAs = index.series.with.NAs
  shifted.series = index.series.with.NAs
  while ( length(index.series.without.NAs[is.na(index.series.without.NAs) == T]) >= 1 ) {
    na.ind = is.na(index.series.without.NAs)
    shifted.series = c(NA, shifted.series[1:(length(shifted.series)-1)])
    index.series.without.NAs[na.ind == T] = shifted.series[na.ind == T]
  }
  index.series.without.NAs
}

# browser()


if ( daily.availability == 1 ) {
  # reduced.temp.data.daily.list = list()
  reduced.temp.data.daily.list[['date']] = reduced.temp.data.daily.list[['date']]
  reduced.temp.data.daily.list[['weekday']] = reduced.temp.data.daily.list[['weekday']]
  for ( i in 1:length(file.names)) {
    reduced.temp.data.daily.list[[file.names[i]]] =
      eliminate.NAs(reduced.temp.data.daily.list[[file.names[i]]])
  }
  if (verbose == T) {
    cat('If daily.availability = 1: NAs have been filled with previous values.\n')
  }
}



# compute daily returns, store everything in a new list:
compute.simple.returns = function(series) {
  ret = diff(series)/series[1:(length(series) - 1)]*100
}

compute.log.returns = function(series) {
  ret = diff(log(series))/log(exp(1))
}

data.daily.list = list()
if ( return.formula == 'simple' ) {
  data.daily.list[['date']] = reduced.temp.data.daily.list[['date']][-1]
  data.daily.list[['weekday']] = reduced.temp.data.daily.list[['weekday']][-1]
  for ( i in 1:length(file.names)) {
    data.daily.list[[file.names[i]]] =
      reduced.temp.data.daily.list[[file.names[i]]][-1]
  }
  # now attach all returns to the list:
  for ( i in 1:length(file.names)) {
    data.daily.list[[paste(file.names[i], '.ret', sep = '')]] =
      compute.simple.returns(reduced.temp.data.daily.list[[file.names[i]]])
  }

}


######################################################################
###  NOW WE PROCEED TO BUILDING THE TWO-DIMENSIONAL SERIES OF WEEKLY VALUES.
######################################################################
if (verbose == T) {cat('Beginning to build the weekly series.\n')}
###  Using temp.data.daily.list, build temp.data.weekly.list successively.
###  First, substitute Wednesdays (or Mondays) for those Tuesdays which are missing.
###  The following two functions do this.
subst.Wednesday = function(weekday, series) {
  new.series = series
  Tuesday.missing = ( is.na(series) & (weekday == '2') ) # Tuesday.missing = ( is.na(series) & (weekday == 'Tuesday') )
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

subst.Monday = function(weekday, series) {
  new.series = series
  Tuesday.missing = ( is.na(series) & (weekday == '2') ) # Tuesday.missing = ( is.na(series) & (weekday == 'Tuesday') )
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


temp.data.weekly.list = list()
temp.data.weekly.list[['date']] = temp.data.daily.list[['date']]
temp.data.weekly.list[['weekday']] = temp.data.daily.list[['weekday']]
for ( i in 1:length(file.names)) {
  temp.subst.series =
    subst.Wednesday(temp.data.daily.list[['weekday']], temp.data.daily.list[[i + 2]])
  temp.subst.series =
    subst.Monday(temp.data.daily.list[['weekday']], temp.subst.series)
  temp.data.weekly.list[[file.names[i]]] = temp.subst.series
}

###  Remove all non-Tuesdays from this list and build it once again
###  (to get the counting index straight):
is.Tuesday = ( temp.data.daily.list[['weekday']] == '2' )  # is.Tuesday = ( temp.data.daily.list[['weekday']] == 'Tuesday' )
temp.data.weekly.list[['date']] = temp.data.weekly.list[['date']][is.Tuesday == T]
temp.data.weekly.list[['weekday']] = temp.data.weekly.list[['weekday']][is.Tuesday == T]
for ( i in 1:length(file.names)) {
  temp.data.weekly.list[[file.names[i]]] = temp.data.weekly.list[[file.names[i]]][is.Tuesday==T]
}


###  depending on weekly.availability, decide which weeks have to be dismissed:
# week.available.1 is TRUE if and only if at least one week's index is available:
week.available.1 = rep(FALSE, length(temp.data.weekly.list[[3]]))
for ( i in 1:length(file.names) ) {
  week.available.1 = (week.available.1 | !is.na(temp.data.weekly.list[[2 + i]]))
}
# week.available.2 is TRUE if and only if all week's indices are available:
week.available.2 = rep(TRUE, length(temp.data.weekly.list[[3]]))
for ( i in 1:length(file.names) ) {
  week.available.2 = (week.available.2 & !is.na(temp.data.weekly.list[[2 + i]]))
}

# In any case, we dismiss all days before the first day for which BOTH indices are available:
# available.from = min(which(day.available.2 == T))
if (weekly.availability == 1) {week.available = week.available.1}
if (weekly.availability == 2) {week.available = week.available.2}
# day.available =
#    c(rep(FALSE, available.from - 1), day.available[available.from:length(day.available)])
if (verbose == T) {cat('Decision which weeks to be dismissed has been made.\n')}

###  Record which weeks are missing:
#   TO BE DONE!!!!!!!!!!!!!!!

###  Now dismiss those weeks for which no data are available:
reduced.temp.data.weekly.list = list()
reduced.temp.data.weekly.list[['date']] = temp.data.weekly.list[['date']][week.available == T]
reduced.temp.data.weekly.list[['weekday']] =
    temp.data.weekly.list[['weekday']][week.available == T]
for ( i in 1:length(file.names)) {
  reduced.temp.data.weekly.list[[file.names[i]]] =
    temp.data.weekly.list[[file.names[i]]][week.available == T]
}
if (verbose == T) {cat('Not available weeks have been dismissed.\n')}

if ( weekly.availability == 1 ) {
  for ( i in 1:length(file.names)) {
    reduced.temp.data.weekly.list[[file.names[i]]] =
      eliminate.NAs(reduced.temp.data.weekly.list[[file.names[i]]])
  }
}


data.weekly.list = list()
if ( return.formula == 'simple' ) {
  data.weekly.list[['date']] = reduced.temp.data.weekly.list[['date']][-1]
  data.weekly.list[['weekday']] = reduced.temp.data.weekly.list[['weekday']][-1]
  for ( i in 1:length(file.names)) {
    data.weekly.list[[file.names[i]]] =
      reduced.temp.data.weekly.list[[file.names[i]]][-1]
  }
  # now attach all returns to the list:
  for ( i in 1:length(file.names)) {
    data.weekly.list[[paste(file.names[i], '.ret', sep = '')]] =
      compute.simple.returns(reduced.temp.data.weekly.list[[file.names[i]]])
  }
}

#####################################################################################
###   Save daily and weekly data in files:
#####################################################################################

file.name.daily = 'data'
for ( i in 1:length(file.names) ) {
  file.name.daily = paste(file.name.daily, '_', file.names[i], sep = '')
}
file.name.daily = paste(file.name.daily, '_daily.dat', sep = '')

file.name.weekly = 'data'
for ( i in 1:length(file.names) ) {
  file.name.weekly = paste(file.name.weekly, '_', file.names[i], sep = '')
}
file.name.weekly = paste(file.name.weekly, '_weekly.dat', sep = '')


# save daily data to file:
output.matrix.daily = c(data.daily.list[[1]], data.daily.list[[2]])
# add level series to the output matrix:
for ( i in 3:(2 + length(file.names)) ) {
  output.matrix.daily = c(output.matrix.daily, data.daily.list[[i]])
}
# add return series to the output matrix:
for ( i in (3 + length(file.names)):(2 + 2*length(file.names)) ) {
  output.matrix.daily = c(output.matrix.daily, data.daily.list[[i]])
}
dim(output.matrix.daily) = c(length(data.daily.list[[3]]), 2+2*length(file.names))

headline = c('date', 'weekday', file.names, paste(file.names, '.ret', sep = ''))
write(headline,  file = file.name.daily, ncolumns = 2+2*length(file.names))
write(t(output.matrix.daily), file = file.name.daily, ncolumns = 2+2*length(file.names),
            append = T)

# save weekly data to file:
output.matrix.weekly = c(data.weekly.list[[1]], data.weekly.list[[2]])
# add level series to the output matrix:
for ( i in 3:(2 + length(file.names)) ) {
  output.matrix.weekly = c(output.matrix.weekly, data.weekly.list[[i]])
}
# add return series to the output matrix:
for ( i in (3 + length(file.names)):(2 + 2*length(file.names)) ) {
  output.matrix.weekly = c(output.matrix.weekly, data.weekly.list[[i]])
}
dim(output.matrix.weekly) = c(length(data.weekly.list[[3]]), 2+2*length(file.names))

headline = c('date', 'weekday', file.names, paste(file.names, '.ret', sep = ''))
write(headline,  file = file.name.weekly, ncolumns = 2+2*length(file.names))
write(t(output.matrix.weekly), file = file.name.weekly, ncolumns = 2+2*length(file.names),
            append = T)

cdts <- list(
  file.names = file.names,
  data.daily = data.daily.list,
  data.weekly = data.weekly.list,
  all.data.daily = temp.data.daily.list)

cat('Total time required: ', difftime(Sys.time(), start.time), '\n')

class(cdts) = "cdts"
cat("Class attributes can be called via the following names:\n")
cat(names(cdts), "\n")
return(cdts)
}
