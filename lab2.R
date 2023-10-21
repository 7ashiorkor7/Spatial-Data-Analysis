#Lab on  Variogram Model fitting, Spatial Regression

#Question3
install.packages("geoR")
library("geoR")
rain <- read.table("RainIowa.txt")
head(rain)
rain.g <- as.geodata(cbind(rain$Latitude, rain$Longitude, rain$mm))

rain.var <- variog(rain.g)
rain.var$max.dist

rain.var <- variog(rain.g, max.dist=3.127459, pairs.min=30)
plot(rain.var)

rain.var1 <- variog(rain.g, estimator.type="modulus", max.dist=3.127459, pairs.min=30)
lines(rain.var1, col=6)

r1 <- lm(mm~Latitude+Longitude, data=rain)
rainfall <- as.geodata(cbind(rain$Latitude, rain$Longitude, r1$residuals))

rain.var2 <- variog(rainfall, max.dist=3.127459, pairs.min=30)
plot(rain.var2)

rain.var3 <- variog(rainfall, estimator.type="modulus", max.dist=3.127459, pairs.min=30)
lines(rain.var3, col=6)

rain.var4 <- variog(rainfall, max.dist=3.12759)

dev.new()
plot(rain.var4)

rain.var <- variog(rainfall, direction=pi/4, max.dist=3.127459)

ini.vals <- expand.grid(seq(0.1, 2000, l=10), seq(0.1,2000,l=10))

rain.wnls <- variofit(rain.var2, ini=ini.vals, fix.nug=FALSE, nugget=1000)
rain.reml <- likfit(rainfall, fix.nug=FALSE, nugget=500, ini.cov.pars=c(1000,1), lik.method="REML", cov.model="matern", kappa=0.5)

plot(rain.var2)
lines(rain.wnls, col=6)
lines(rain.reml, col=5)


#Question4
soil <- read.table("soilph.txt")
colnames(soil) <- c("x","y","ph")
r <- soil$x
c <- soil$y
r2 <- soil$x*soil$y
rc <- soil$x*soil$y
c2 <- soil$y*soil$y
soil.t <- lm(ph~r+c+r2+rc+c2, data=soil)
resid.data <-as.data.frame(cbind(soil$x,soil$y,resid=soil.t$resid))
soil.res <- as.geodata(resid.data)

soil.var <- variog(soil.res)
soil.var <- variog(soil.res, max.dist=7.071068, min.pair=30)

soil.robust <- variog(soil.res, estimator.type="modulus", max.dist=7.071068, min.pair=30)
plot(soil.var, ylim=c(0,0.045))
points(soil.robust$u,soil.robust$v, col=2)

soil.dir <- variog4(soil.res, max.dist=7.071068, min.pair=30)
plot(soil.dir)

variog_N <- variog(soil.res, max.dist =7.071068, direction=0, tolerance = 15, unit.angle="degrees")
variog_NE <- variog(soil.res, max.dist =7.071068, direction=45, tolerance = 15, unit.angle="degrees")
variog_E <- variog(soil.res, max.dist =7.071068, direction=90, tolerance = 45, unit.angle="degrees")
variog_SE <- variog(soil.res, max.dist =7.071068, direction=135, tolerance = 45, unit.angle="degrees")

dev.new()
plot(variog_N, type="b", xlim=c(0,8), ylim=c(0,0.05))
lines(variog_NE, type="b", col="red")
lines(variog_E, type="b", col="blue")
lines(variog_SE, type="b", col="green")

ini.vals <- expand.grid(seq(0.005, 0.1,l=10), seq(0.01,5,l=10))
soil.ols <- variofit(soil.var, ini=ini.vals,fix.nug=FALSE, max.dist=7.071068, wei="equal")
soil.wnls <- variofit(soil.var, ini=ini.vals,fix.nug=FALSE, wei="cressie", max.dist=7.071068)
plot(soil.var)
lines(soil.ols, col=4)
lines(soil.wnls, col=6, lty=2)

soilph <- as.geodata(soil)
soil.ml <- likfit(soilph, trend="2nd", ini=c(1.5,1.2), fix.nug=FALSE)
soil.reml <- likfit(soilph, trend="2nd", ini=c(1.5,1.2), fix.nug=FALSE)

plot(soil.var, ylim=c(0,0.05))
lines(soil.ml, lty=2)

ols.dat <- as.data.frame(cbind(coords=soilph$coords, resid=soil.t$resid))
soil.ols <- as.geodata(ols.dat)
X <- cbind(rep(1,121), r,c,r2,c2)
soil.var1 <- variog(soil.ols, max.dist = 7.071068)
soil.wnls1 <- variofit(soil.var1, ini=ini.vals, fix.nug=TRUE, fit.nig=TRUE, max.dist = 7.071068)

dist.mat <- matrix(scan("distance.txt"), ncol=121, byrow = T)
dist.mat[1:4,1:4]

Vhat <- matrix(NA, nrow=121, ncol=121)
Vhat <- cov.spatial(dist.mat, cov.model="exp", cov.pars = soil.wnls1$cov.pars)
beta.hat1 <- solve(t(X)%*%solve(Vhat)%*%X)%*%t(X)%*%solve(Vhat)%*%soil$ph

