library(ggplot2)
library(reshape)
library(zoo)

dat <- read.csv("~/Documents/Projects_Current/partitionfinder/docs/downloadstats.csv")
dat$date <- as.yearmon(dat$date)
dat <- aggregate(. ~ date, data=dat, FUN=sum, na.action=na.pass)
dat <- dat[order(as.Date(dat$date)),]

dat$Mac <- cumsum(dat$Mac)
dat$Windows <- cumsum(dat$Windows)
 

# Download plot
m <- dat[,c(1,2)]
m$OS <- "Mac"
names(m) <- c('date', 'downloads', 'OS')
w <- dat[,c(1,3)]
w$OS <- "Windows"
names(w) <- c('date', 'downloads', 'OS')

dl <- rbind(m, w)
dl$date <- as.Date(dl$date, format = "%d/%m/%Y")

quartz(width=10, height=5)
p <- ggplot(dl, aes(date, downloads))
p + geom_area(aes(colour = OS, fill= OS), position = 'stack')
dev.copy2pdf(file="~/Documents/Projects_Current/partitionfinder/docs/downloadplot.pdf")
dev.off()


# Citation plot
c <- dat[,c(1,5)]
c$date <- as.Date(c$date, format = "%d/%m/%Y")
c <- c[complete.cases(c),]

quartz(width=10, height=5)
p <- ggplot(c, aes(date, citations.googlescholar))
p + geom_point() + geom_smooth(stat="identity")
dev.copy2pdf(file="~/Documents/Projects_Current/partitionfinder/docs/citationsplot.pdf")
dev.off()

