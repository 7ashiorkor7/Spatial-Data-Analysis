#Lab on Kriging

install.packages("geoR")
library("geoR")

#Question2
#a
c0 <- 0
c1 <- 1
phi <- 4

c0 <- 0.25
c1 <- 0.75
phi <- 2

xc <- c(0,1,2,3,4,3,2)
yc <- c(0,1,3,2,1,0,1)

dd <- cbind(xc,yc)
m <- length(xc)

dist <- matrix(0, nrow = m, ncol = m)
for (i in 1:m) {
  for (j in 1:m) {
    dist[i,j] = ((dd[i,1]-dd[j,1])^2 + (dd[1,2]-dd[j,2]^2)^.5)
    
  }
  
}

#Computing the Gamma matrix for kriging system of equations
G <- matrix(NA, nrow=m, ncol=m)
for(i in 1:m-1){
  for(j in 1:m-1){
    G[i,j] = c1*exp(-dist[i,j]/phi)
  }
}

G[m,] <- rep(1,m)
G[,m] <- rep(1,m)
G[m,m] <- 0

D <- rep(0,m)
  for(j in 1:m){
    D[j] = c1*exp(-dist[m,j]/phi)
  }
  D[m] = 1
  
wt <- solve(G)%*%D
wt <- round(solve(G) %*% D, digits=3)
wt

var_z_hat <- round(c0+c1-t(wt) %*% D, digits=3)
var_z_hat

plot(xc[-7], yc[-7], xlab="X coordinate", ylab="Y coordinate", main="Kriging weights", pch=16) 

points(xc[7], yc[7], pch=13)

text(0,0, label= wt[1], pos=4) #plot the weights at each location
text(1,1, label= wt[2], pos=3) 
text(2,3, label= wt[3], pos=1) 
text(3,2, label= wt[4], pos=3)  
text(4,1, label= wt[5], pos=2) 
text( 3.0, label= wt[6], pos=3) 

#b  
c0 <- 0
c1 <- 1
phi <- 4

c0 <- 0.25
c1 <- 0.75
phi <- 2

x=xc[1:6]
y=yc[1:6]
n <- length(x)
m = n+3

G <- matrix(NA, nrow = m, ncol = m)  
  for(i in 1:n){
    for(j in 1:n){
      G[i,j] = c1*exp(-dist[i,j]/phi)
    }
  }

G[m,] <- rep(1,m)
G[,m] <- rep(1,m)
G[m,m] <- 0

G[1:6, n+1]=x
G[1:6, n+2]=y
G[n+1, 1:6]=x
G[n+2, 1:6]=y
G[,m] <- rep(1,m)
G[m,] <- rep(1,m)
diag(G)=0
G[7:9, 7:9]=0

g <- rep(0,m)
  for(j in 1:n){
    g[j] = c1 * exp(-dist[n+1,j]/phi)
  }
  g[n+1] = 2
  g[n+2] = 1
  g[m] = 1
  
wt <- solve(G) %*% g
wt <- round(solve(G) %*%g, digits = 3)
wt

################################################
  
#Question3

soilph.d <- read.table('soilph.txt')
soilph <- as.geodata(soilph.d)
  
newdata1 <- c(1.5, 1.5)

soil.kr1 <- krige.conv(soilph, locations=newdata1, krige=krige.control(type.krige="ok",trend.d="2nd",trend.l="2nd",cov.pars=c(0.0363,0.6143),nugget=0.009,cov.model="exponential"))

soil.kr1$predict
soil.kr1$krige.var

soil.ok1 <- krige.conv(soilph, locations=newdata1, krige=krige.control(type.krige="ok",cov.pars=c(0.0363,0.6143), nugget=0.009, cov.model="exponential"))

soil.ok1$predict
soil.ok1$krige.var

newdata2<-c(5.5, 5.5)
soil.kr2<-krige.conv(soilph, locations=newdata2, krige=krige.control(type.krige="ok",trend.d="2nd", trend.l="2nd", cov.pars=c(0.0333,0.6143),nugget=0.009, cov.model="exponential"))

soil.kr2$predict
soil.kr2$krige.var

newdata3 <- c(8,3.25)
soil.kr3<-krige.conv(soilph, locations=newdata3, krige=krige.control(type.krige="ok",
trend.d="2nd", trend.l="2nd", cov.pars=c(0.0333,0.6143),nugget=0.009, cov.model="exponential"))

soil.kr3$predict
soil.kr3$krige.var

soil.reml <- likfit(soilph,trend="2nd", ini=c(0.3,1), fix.nug=TRUE, met="REML")
soil.reml
scross.reml <- xvalid(soilph, model=soil.reml)

soil.var <- variog(soilph,trend="2nd", max.dist=7.071068,min.pair=30)
soil.wnls <- variofit(soil.var, ini=c(0.05,1), fix.nug=TRUE) 
soil.cross.wnls <- xvalid(soilph, model=soil.wnls)

mean(scross.reml$error^2)
mean(scross.reml$std)
sqrt(mean(scross.reml$std^2))

mean(soil.cross.wnls$error^2)
mean(soil.cross.wnls$std)

