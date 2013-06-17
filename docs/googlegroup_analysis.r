d <- read.csv("~/Documents/Projects_Current/partitionfinder/docs/googlegroup_stats.csv")

format = "%a, %d %b %Y %H:%M:%S"

d$Sent <- strptime(d$Sent, format=format)
d$Replied <- strptime(d$Replied, format=format)

d$ResponseTime <- d$Replied - d$Sent

median(d$ResponseTime)
