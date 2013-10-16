library(ggplot2)
library(reshape)

dat <- read.csv("~/Documents/Projects_Current/partitionfinder/docs/downloadstats.csv")
dat <- aggregate(. ~ date, data=dat, FUN=sum, na.action=na.pass)
dat <- dat[order(as.Date(dat$date, format="%d/%m/%Y")),]

dat$Mac <- cumsum(dat$Mac)
dat$Windows <- cumsum(dat$Windows)
 

# Download plot
dl <- melt(dat[,c(1, 2, 3)])
colnames(dl) <- c("date", "OS", "downloads")
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

