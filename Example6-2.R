####################################
#  Code for Simulations in Chapter 12.6.2
####################################

## Load algorithm functions
source('HelperFunctions.R')
library(rgl)

##################################################################
## First scenario, in which spectral clustering performs very well
##################################################################

## Simulate Data from a linear factor model with noise

set.seed(100)
d <- 3
b1 <- c(0.1,0.1,0.8)
b2 <- c(0.1,0.8,0.1)
b3 <- c(0.8,0.1,0.1)

a1 <- b1/sqrt(sum(b1^2))
a2 <- b2/sqrt(sum(b2^2))
a3 <- b3/sqrt(sum(b3^2))

n <- 10000
d <- 3
q <- 0.95

Z1 <- (-log(runif(n)))^(-1)
Z2 <- (-log(runif(n)))^(-1)
Z3 <- (-log(runif(n)))^(-1)
eta <- (-log(runif(n)))^(-1)
N1 <- abs(rnorm(n))
N2 <- abs(rnorm(n))
N3 <- abs(rnorm(n))
X <- matrix(NA,n,d)
for (i in 1:n){
  X[i,] <- b1*Z1[i]+ b2*Z2[i] + b3*Z3[i] + 0.5*eta[i] * c(N1[i],N2[i],N3[i])
}

# Process data

# Transform and truncate
Frechettrans<-function(x) 1/(1-ecdf(x)(x)*length(x)/(length(x)+1))
X.ext <- apply(X,2,Frechettrans)
norm_vec <- function(x) sqrt(sum(x^2))
norms<-apply(X.ext,1,norm_vec)

X.ext <- X.ext[norms>quantile(norms,q),]
ang.ext <- X.ext/norms[norms>quantile(norms,q)]


###########################
# k-means
###########################

library(skmeans)
k <- 3
fit.kmeans <- clusterMeans(ang.ext,k)

# Plot results

view<-cbind(c(1,0,0,0),c(0,0.33,-0.93,0),c(0,0.93,0.33,0),c(0,0,0,1))
open3d(windowRect = 0 + c( 0, 0, 600,600 ) )
par3d(userMatrix = view)
plot3d(ang.ext[fit.kmeans$assignments==1,],
       xlab = "", ylab = "", zlab = "",
       axes = FALSE,col='red',cex=1.9,pch=0)
plot3d(ang.ext[fit.kmeans$assignments==2,],cex=1.9,col="blue",pch=1,add=T)
plot3d(ang.ext[fit.kmeans$assignments==3,],cex=1.9,pch=1,add=T)
arc3d(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(1, 0, 0), c(0, 0, 1), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(0, 0, 1), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
aspect3d(1, 1, 1)
axes3d(edges = c("x-+", "y-+", "z+-"),
       ntick = 6,                       # Attempt 6 tick marks on each side
       cex = 1)                       # Smaller font
# Indicate points of mass
scale<-seq(0:20)/20
plot3d(scale%*%t(a1),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a2),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a3),cex=1,pch=1,add=T,col="darkgray")
# Fix angle
pp <- dget("perspective.R")
par3d(pp)

# Save
rgl.snapshot("LinNoisesphericalV1.png", fmt="png")


###########################
# k-pc
###########################
# source('KpcFns.R')
centers.kpc <- clusterPC(ang.ext,k)
centers.kpc
assignments.kpc <- getClusterIndex(X.ext,centers.kpc)

open3d(windowRect = 0 + c( 0, 0, 600,600 ) )
par3d(userMatrix = view)
plot3d(ang.ext[assignments.kpc==1,],
       xlab = "", ylab = "", zlab = "",
       axes = FALSE,col='red',cex=1.9,pch=0)
plot3d(ang.ext[assignments.kpc==2,],cex=1.9,pch=1,add=T)
plot3d(ang.ext[assignments.kpc==3,],col="blue",cex=1.9,pch=1,add=T)
arc3d(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(1, 0, 0), c(0, 0, 1), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(0, 0, 1), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
aspect3d(1, 1, 1)
axes3d(edges = c("x-+", "y-+", "z+-"),
       ntick = 5,                       # Attempt 6 tick marks on each side
       cex = 0.9)                       # Smaller font

# Indicate points of mass
plot3d(scale%*%t(a1),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a2),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a3),cex=1,pch=1,add=T,col="darkgray")

# Fix angle
par3d(pp)

# This leads to exactly the same graphic as the spherical K-means,
# so we do not save this plot

###########################
# spectral clustering
###########################

# Version 1: k=3, i.e. number of clusters equals number of factors

assignments.spectral <- clusterSpectral(ang.ext,k,10)


open3d(windowRect = 0 + c( 0, 0, 600,600 ) )
par3d(userMatrix = view)
plot3d(ang.ext[assignments.spectral==2,],
       xlab = "", ylab = "", zlab = "",
       axes = FALSE,col='blue',cex=1.9,pch=0)
plot3d(ang.ext[assignments.spectral==1,],cex=1.9,pch=1,add=T)
plot3d(ang.ext[assignments.spectral==3,],col="red",cex=1.9,pch=1,add=T)
arc3d(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(1, 0, 0), c(0, 0, 1), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(0, 0, 1), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
aspect3d(1, 1, 1)
axes3d(edges = c("x-+", "y-+", "z+-"),
       ntick = 6,                       # Attempt 6 tick marks on each side
       cex = 1) # Smaller font

# Indicate points of mass
plot3d(scale%*%t(a1),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a2),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a3),cex=1,pch=1,add=T,col="darkgray")

# Fix angle
par3d(pp)

# Save
rgl.snapshot("LinNoiseSpectralk3V1.png", fmt="png")


# Version 2: k=4, i.e. number of clusters is one larger than number of factors

assignments.spectral <- clusterSpectral(ang.ext,k+1,10)

open3d(windowRect = 0 + c( 0, 0, 600,600 ) )
plot3d(ang.ext[assignments.spectral==1,],
       xlab = "", ylab = "", zlab = "",
       axes = FALSE,col='chartreuse4',cex=1.9)
plot3d(ang.ext[assignments.spectral==2,],cex=1.9,col="blue",pch=1,add=T)
plot3d(ang.ext[assignments.spectral==3,],cex=1.9, pch=1,add=T)
plot3d(ang.ext[assignments.spectral==4,],cex=1.9,col="red", pch=1,add=T)
arc3d(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(1, 0, 0), c(0, 0, 1), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(0, 0, 1), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
aspect3d(1, 1, 1)
axes3d(edges = c("x-+", "y-+", "z+-"),
       ntick = 6,                       # Attempt 6 tick marks on each side
       cex = 1)                       # Smaller font

# Indicate points of mass
plot3d(scale%*%t(a1),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a2),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a3),cex=1,pch=1,add=T,col="darkgray")

# Fix angle
par3d(pp)

# Save
rgl.snapshot("LinNoiseSpectralk4V1.png", fmt="png")

##############################################################################
## Second scenario, in which spectral clustering with k=4 performs not so well
##############################################################################

## Simulate Data from a linear factor model with noise

set.seed(100)
d <- 3
b1 <- c(0.25,0.5,0.25)
b2 <- c(0.25,0.25,0.5)
b3 <- c(0.5,0.25,0.25)
a1 <- b1/sqrt(sum(b1^2))
a2 <- b2/sqrt(sum(b2^2))
a3 <- b3/sqrt(sum(b3^2))

n <- 10000
d <- 3
q <- 0.95

Z1 <- (-log(runif(n)))^(-1)
Z2 <- (-log(runif(n)))^(-1)
Z3 <- (-log(runif(n)))^(-1)
eta <- (-log(runif(n)))^(-1)
N1 <- abs(rnorm(n))
N2 <- abs(rnorm(n))
N3 <- abs(rnorm(n))
X <- matrix(NA,n,d)
for (i in 1:n){
  X[i,] <- b1*Z1[i]+ b2*Z2[i] + b3*Z3[i] + 0.5*eta[i] * c(N1[i],N2[i],N3[i])
}

# Process data

# Transform and truncate
Frechettrans<-function(x) 1/(1-ecdf(x)(x)*length(x)/(length(x)+1))
X.ext <- apply(X,2,Frechettrans)
norm_vec <- function(x) sqrt(sum(x^2))
norms<-apply(X.ext,1,norm_vec)

X.ext <- X.ext[norms>quantile(norms,q),]
ang.ext <- X.ext/norms[norms>quantile(norms,q)]


###########################
# k-means
###########################

library(skmeans)
k <- 3
fit.kmeans <- clusterMeans(ang.ext,k)

# Plot results

view<-cbind(c(1,0,0,0),c(0,0.33,-0.93,0),c(0,0.93,0.33,0),c(0,0,0,1))
open3d(windowRect = 0 + c( 0, 0, 600,600 ) )
par3d(userMatrix = view)
plot3d(ang.ext[fit.kmeans$assignments==1,],
       xlab = "", ylab = "", zlab = "",
       axes = FALSE,col='red',cex=1.9,pch=0)
plot3d(ang.ext[fit.kmeans$assignments==2,],cex=1.9,col="blue",pch=1,add=T)
plot3d(ang.ext[fit.kmeans$assignments==3,],cex=1.9,pch=1,add=T)
arc3d(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(1, 0, 0), c(0, 0, 1), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(0, 0, 1), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
aspect3d(1, 1, 1)
axes3d(edges = c("x-+", "y-+", "z+-"),
       ntick = 6,                       # Attempt 6 tick marks on each side
       cex = 1)                       # Smaller font
# Indicate points of mass
scale<-seq(0:20)/20
plot3d(scale%*%t(a1),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a2),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a3),cex=1,pch=1,add=T,col="darkgray")
# Fix angle
pp <- dget("perspective.R")
par3d(pp)

# Save
rgl.snapshot("LinNoisesphericalV2.png", fmt="png")


###########################
# k-pc
###########################
# source('KpcFns.R')
centers.kpc <- clusterPC(ang.ext,k)
centers.kpc
assignments.kpc <- getClusterIndex(X.ext,centers.kpc)

open3d(windowRect = 0 + c( 0, 0, 600,600 ) )
par3d(userMatrix = view)
plot3d(ang.ext[assignments.kpc==1,],
       xlab = "", ylab = "", zlab = "",
       axes = FALSE,col='blue',cex=1.9,pch=0)
plot3d(ang.ext[assignments.kpc==2,],col="red",cex=1.9,pch=1,add=T)
plot3d(ang.ext[assignments.kpc==3,],,cex=1.9,pch=1,add=T)
arc3d(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(1, 0, 0), c(0, 0, 1), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(0, 0, 1), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
aspect3d(1, 1, 1)
axes3d(edges = c("x-+", "y-+", "z+-"),
       ntick = 5,                       # Attempt 6 tick marks on each side
       cex = 0.9)                       # Smaller font

# Indicate points of mass
plot3d(scale%*%t(a1),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a2),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a3),cex=1,pch=1,add=T,col="darkgray")

# Fix angle
par3d(pp)

# This leads to exactly the same graphic as the spherical K-means,
# so we do not save this plot

###########################
# spectral clustering
###########################

# Version 1: k=3, i.e. number of clusters equals number of factors

assignments.spectral <- clusterSpectral(ang.ext,k,10)


open3d(windowRect = 0 + c( 0, 0, 600,600 ) )
par3d(userMatrix = view)
plot3d(ang.ext[assignments.spectral==2,],
       xlab = "", ylab = "", zlab = "",
       axes = FALSE,cex=1.9,pch=0)
plot3d(ang.ext[assignments.spectral==1,],cex=1.9,col='blue',pch=1,add=T)
plot3d(ang.ext[assignments.spectral==3,],col="red",cex=1.9,pch=1,add=T)
arc3d(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(1, 0, 0), c(0, 0, 1), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(0, 0, 1), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
aspect3d(1, 1, 1)
axes3d(edges = c("x-+", "y-+", "z+-"),
       ntick = 6,                       # Attempt 6 tick marks on each side
       cex = 1) # Smaller font

# Indicate points of mass
plot3d(scale%*%t(a1),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a2),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a3),cex=1,pch=1,add=T,col="darkgray")

# Fix angle
par3d(pp)

# Save
rgl.snapshot("LinNoiseSpectralk3V2.png", fmt="png")


# Version 2: k=4, i.e. number of clusters is one larger than number of factors

assignments.spectral <- clusterSpectral(ang.ext,k+1,10)

open3d(windowRect = 0 + c( 0, 0, 600,600 ) )
plot3d(ang.ext[assignments.spectral==1,],
       xlab = "", ylab = "", zlab = "",
       axes = FALSE,col='red',cex=1.9)
plot3d(ang.ext[assignments.spectral==2,],cex=1.9,pch=1,add=T)
plot3d(ang.ext[assignments.spectral==3,],cex=1.9,col="blue", pch=1,add=T)
plot3d(ang.ext[assignments.spectral==4,],cex=1.9,col="chartreuse4", pch=1,add=T)
arc3d(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(1, 0, 0), c(0, 0, 1), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(0, 0, 1), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
aspect3d(1, 1, 1)
axes3d(edges = c("x-+", "y-+", "z+-"),
       ntick = 6,                       # Attempt 6 tick marks on each side
       cex = 1)                       # Smaller font

# Indicate points of mass
plot3d(scale%*%t(a1),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a2),cex=1,pch=1,add=T,col="darkgray")
plot3d(scale%*%t(a3),cex=1,pch=1,add=T,col="darkgray")

# Fix angle
par3d(pp)

# Save
rgl.snapshot("LinNoiseSpectralk4V2.png", fmt="png")