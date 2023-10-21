#Lab on Visualization, exploratory analysis of geostatistical data
install.packages("maps") 
library(maps)

rain <- read.table("RainIowa.txt")
par(mfrow=c(1,1))

plot(rain$Lon, rain$Lat, xlim=c(-97,-90), ylim=c(40,44), xlab="Longitude",
     ylab="Latitude", main="Rain Locations in Iowa", "n")
map("county", "iowa", add=TRUE)
points(rain$Lon, rain$Lat, cex=rain$mm/400)

rain_col <- c("yellow", "orange", "red", "green", "black")
rain_levels2 <- cut(rain$mm, c(700, 770, 789.25, 908.50, 975.75, 1043.00))
as.numeric(rain_levels2)

plot(rain$Lon, rain$Lat, xlim=c(-97,-90), ylim=c(40,44), xlab="Longitude",
     ylab="Latitude", main="Rain locations in Iowa", "n")
map("county", "iowa", add=TRUE)
points(rain$Lon, rain$Lat, cex=rain$mm/mean(rain$mm), col=rain_col[as.numeric(rain_levels2)], pch=19)

