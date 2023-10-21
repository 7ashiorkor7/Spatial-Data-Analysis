#Lab on Areal data analysis with Moran's I, Geary's c and autoregressive models

#Question1
install.packages("sp")
install.packages("spdep")
install.packages("spatialreg")
install.packages("maptools")

library(sf)
library(geoR)
library(sp)
library(gstat)
library(spdep)
library(spatialreg)
library(maptools)

acex1 <- read.table('acex1.txt')
acex2 <- read.table('acex2.txt')
acex3 <- read.table('acex3.txt')

temp <- expand.grid(x=1:5, y=1:5)
acex1.sp <- data.frame(temp, z = acex1$V2)
acex2.sp <- data.frame(temp, z = acex2$V2)
acex3.sp <- data.frame(temp, z = acex3$V2)

coordinates(acex1.sp) <- c('x','y')
gridded(acex1.sp) <- T

coordinates(acex2.sp) <- c('x','y')
gridded(acex2.sp) <- T

coordinates(acex3.sp) <- c('x','y')
gridded(acex3.sp) <- T

par(mfrow=c(1,3))
spplot(acex1.sp)
spplot(acex2.sp)
spplot(acex3.sp)

acex_kn1 <- knn2nb(knearneigh(coordinates(acex1.sp), k=1))
acex_kn1_w <- nb2listw(acex_kn1, style="W")
acex_kn1_b <- nb2listw(acex_kn1, style="B")
acex.nb <- cell2nb(5,5,type="rook")
acex_b <- nb2listw(acex.nb, style="B")

acex1.mo <- moran.test(acex1$V2, acex_kn1_w, alternative="two.sided")
acex1.gr <- geary.test(acex1$V2, listw=acex_kn1_w)

acex2.mo <- moran.test(acex2$V2, acex_kn1_w, alternative="two.sided")
acex2.gr <- geary.test(acex2$V2, listw=acex_kn1_w)


acex3.mo <- moran.test(acex3$V2, acex_kn1_w, alternative="two.sided")
acex3.gr <- geary.test(acex3$V2, listw=acex_kn1_w)


###################################################################################################

#Question2
load("scotland.Rdata")
head(scotland)

spplot(scotland, "logratio")
spplot(scotland, "logratio", col.regions=terrain.colors(16))

scotland.nb <- poly2nb(scotland)
par(mar=c(3,3,0,0)+0.2, mgp=c(2,0.8,0))
plot(scotland)
plot(scotland.nb, coordinates(scotland), add=T, col=4)

nneigh <- sapply(scotland.nb, length)
table(nneigh)

scotland.w <- nb2listw(scotland.nb, style = 'B')
scotland.mo <- moran.test(scotland$logratio, scotland.w)

scotland.gr <- geary.test(scotland$logratio, scotland.w)

scotland.mc <- moran.mc(scotland$logratio, scotland.w, 999)

scotland_lm<-lm(percentAFF~logratio,data=scotland)   
summary(scotland_lm)
moran.test(scotland_lm$residuals, listw=scotland.w)

scotland.sar <- errorsarlm( percentAFF ~ 1, data=scotland, scotland.w, method="eigen", quiet=FALSE)
summary(scotland.sar)

scotlandpop.sar <-errorsarlm(percentAFF~logratio, data=scotland, scotland.w, method="eigen", quiet=FALSE)
summary(scotlandpop.sar) 

moran.test(scotlandpop.sar$residuals, listw=scotland.w)

