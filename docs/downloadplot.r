library(ggplot2)
library(reshape)

dat <- read.csv("~/Documents/Projects_Current/partitionfinder/docs/downloadstats.csv")

#the version changes in march and may need to be fixed up
dat <- aggregate(. ~ date, data=dat, FUN=sum)

dat <- dat[order(as.Date(dat$date, format="%d/%m/%Y")),]

dat <- dat[,-c(4)]

dat$Mac <- cumsum(dat$Mac)
dat$Windows <- cumsum(dat$Windows)
 
dat <- melt(dat)

colnames(dat) <- c("date", "OS", "downloads")

dat$date <- as.Date(dat$date, format = "%d/%m/%Y")


p <- ggplot(dat, aes(date, downloads))

p + geom_area(aes(colour = OS, fill= OS), position = 'stack')

