#Lab on simulation of random fields, empirical semivariograms, covariance models

#Question3
soil <- read.table("soilph.txt")
colnames(soil) <- c("x","y","ph")

#a
library(scatterplot3d)
scatterplot3d(soil$x, soil$y, soil$ph)

summary(soil$ph)
sqrt(var(soil[,3]))
stem(soil[,3])
boxplot(soil[,3])

par(mfrow=c(2,2))
soil.mat <- tapply(soil[,3], list(factor(soil[,1]),factor(soil[,2])),function(x)x)
  contour(soil.mat)
  image(soil.mat)
  persp(soil.mat)

plot(soil$x, soil$y, type="n", xlab="n",ylab=" ")
text(soil$x, soil$y, round(soil$ph,1), col=3)

library(geoR)

soil.g <- as.geodata(soil)
plot(soil.g)
points(soil.g)

#b
r <- soil$x
c <- soil$y
r2 <- soil$x*soil$y
rc <- soil$x*soil$y
c2 <- soil$y*soil$y

fit.s <- lm(ph ~ r+c+r2+rc+c2, data=soil)
fit.r <- trend.spatial(soil.g, trend = "2nd")

par(pty="s", mfrow=c(1,2))
soilgrid.mat <- tapply(fit.s$resid, list(factor(soil$x), factor(soil$y)), function(x)x)

plot(soilgrid.mat[,-1], soilgrid.mat[,-11], xlab="column Z", ylab="column Z-1")
abline(c(0,1))
plot(soilgrid.mat[-1,], soilgrid.mat[-11,], xlab="pH in row Z-1")
abline(c(0,1))

#c
sloe <- loess(soil$ph ~ soil$x*soil$y, data=soil, normalize=F,span=.15)
sloe.pred <- predict(sloe)

soil.grid <- as.data.frame(cbind(x=soil$x, y=soil$y))
sloe.mat<-tapply(sloe.pred, list(factor(soil.grid$x), factor(soil.grid$y)), function(x)x)
sloeRes.mat<-tapply(sloe$res,list(factor(soil.grid$x), factor(soil.grid$y)),function(x)x)

plot(sloeRes.mat[,-1], sloeRes.mat[,-11])
abline(c(0,1))
plot(sloeRes.mat[-1,], sloeRes.mat[-11,])
abline(c(0,1))
mtext('h scatter plot of LOESS residuals', line=2)

#d
par(mfrow=c(2,2))
soilquad.mat <- tapply(fit.s$fit, list(factor(soil$x), factor(soil$y)), function(x)x)

image(soilquad.mat)
persp(soilquad.mat)
mtext("Surface fitted by Quadratic function")

image(sloe.mat)
persp(sloe.mat)
mtext("Surface fitted by Quadratic function")

#e
install.packages("gstat")
library(gstat)
soil.mp <- medpolish(soil.mat, na.rm=T)
soil.trend <- soil.mat - soil.mp$residuals

par(mfrow=c(1,3))
image(soil.mat)
image(soil.trend)
image(soil.mp$res)

persp(soil.mat)
persp(soil.trend)
persp(soil.mp$res)

##############################################
#question4
xy16 <- expand.grid(x=seq(1,16,len=16), y=seq(1,16), len=16)

z.s1<-grf(16*16, grid=xy16, cov.model="gaussian", cov.pars = c(1,1))
z.s2<-grf(16*16, grid=xy16, cov.model="gaussian", cov.pars = c(1,3))
z.s3<-grf(16*16, grid=xy16, cov.model="gaussian", cov.pars = c(1,5))
z.s4<-grf(16*16, grid=xy16, cov.model="gaussian", cov.pars = c(1,1), nugget = 6)

ez1.var <- variog(z.s1)
sqrt(2*16^2)/2

ez1.var <- variog(z.s1, max.dist=11.31, min.pair=30)
ez2.var <- variog(z.s2, max.dist=11.31, min.pair=30)
ez3.var <- variog(z.s3, max.dist=11.31, min.pair=30)
ez4.var <- variog(z.s4, max.dist=11.31, min.pair=30)

par(mfrow=c(2,2))
plot(ez1.var)
plot(ez2.var)
plot(ez3.var)
plot(ez4.var)

z.sp1 <- grf(16*16, grid=xy16, cov.model ="spherical", cov.pars = c(1,3))
z.sp2 <- grf(16*16, grid=xy16, cov.model ="spherical", cov.pars = c(1,5))

zs1.v <- variog(z.sp1, max.dist=11.31, min.pair=30)
zs2.v <- variog(z.sp2, max.dist=11.31, min.pair=30)

plot(zs1.v)
plot(zs1.v)

persp(z.s1)
persp(z.s2)

#################################################################################
#question5
s1 = (0,1) s2=(0,0)

x <- c(0,0,1,1)
y <- c(1,0,0,1)

q <- as.data.frame(cbind(x,y))

x1 <- rep(rep(0,4),4)
dist <- matrix(x1, nrow=4, ncol=4)

C <- matrix(x1, nrow=4, ncol=4)

for(i in 1:4){
  for(j in 1:4){
    dist[i,j]=((q[i,1]-q[j,1]^2)+(q[i,2])^2)^.5
    C[i,j]=2.5*exp(-dist[i,j]/2)
  }
}  

a<-c(0.5, 0.2, 0.2, 0.1)
var_z_hat <- t(a) %*% C %*% a
var_z_hat
  

