library(ggplot2)
library(reshape)

dat <- read.csv("~/Documents/Projects_Current/partitionfinder/docs/downloadstats.csv")

dat$Mac <- cumsum(dat$Mac)
dat$Windows <- cumsum(dat$Windows)
 
dat <- melt(dat)

colnames(dat) <- c("date", "ver", "OS", "downloads")

dat$date <- as.Date(dat$date, format = "%d/%m/%Y")


p <- ggplot(dat, aes(date, downloads))

p + geom_area(aes(colour = OS, fill= OS), position = 'stack')

