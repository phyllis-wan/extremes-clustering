#############################################
# Code for Figures in Chapter 12.1
#############################################

library(graphicalExtremes)
library(rgl)
library(plot3D)

#############################################
# Function to standardize margins to unit Pareto
#############################################

Paretotrans<-function(x) 1/(1-ecdf(x)(x)*length(x)/(length(x)+1)) 

#############################################
# 3D plot of river extremes, Plot 1.1
#############################################

# Load danube river data and select Stations 27,28,29
data <- danube$data_clustered
n <- nrow(data)
i1 <- 27
i2 <- 28
i3 <- 29
dat <- data[,c(i1,i2,i3)]

# Transform the data to Pareto margin
river.ext <- apply(dat,2,Paretotrans)

# 3D-plot of transformed sample
scatter3D(river.ext[,1], river.ext[,2], river.ext[,3],colvar = NULL, bty = "b2",col = "black", phi = 5, theta=35, pch = 19, cex = 1.5, xlab = "Station 27",
          ylab ="Station 28", zlab = "Station 29")


#############################################
# Example of linear factor model, Plot 1.2
#############################################

# Generate sample (n=1000) from three-dimensional linear factor model
set.seed(4)
n<-1000
Z1<--1/log(runif(n))
Z2<--1/log(runif(n))
Z3<--1/log(runif(n))
X1<-rep(0,n)
X2<-rep(0,n)
X3<-rep(0,n)
for (i in 1:n){
  X1[i]<-0.07*Z1[i]+0.01*Z2[i]+0.07*Z3[i]
  X2[i]<-0.07*Z1[i]+0.06*Z2[i]+0.01*Z3[i]
  X3[i]<-0.07*Z1[i]+0.07*Z2[i]+0.01*Z3[i]
}
X<-cbind(X1,X2,X3)

# Transform margins for easier comparison

dim(X)<-c(1000,3)
X<-apply(X,2,Paretotrans)

# 3D-plot of transformed sample

scatter3D(X[,1], X[,2], X[,3],colvar = NULL, bty = "b2",col = "black", phi = 5, theta=35, pch = 19, cex = 1.5,xlab = "",
          ylab ="", zlab = "")
text3D(500, 0, -150, labels = expression(X[1]), add = TRUE)
text3D(1050, 600, -150, labels = expression(X[2]), add = TRUE)
text3D(-300, 0, 500, labels = expression(X[3]), add = TRUE)
#############################################
# Spectral plots
#############################################

# Select the extreme observations above a quantile
q <- 0.95
norm_vec <- function(x) sqrt(sum(x^2)) 

# Apply to River data
norms <- apply(river.ext,1,norm_vec)
river.ext <- river.ext[norms>quantile(norms,q),]
norms <- apply(river.ext,1,norm_vec)
river.ext <- river.ext/norms 

# Plot 3D plot of the extreme observation
open3d(windowRect = 20 + c(0,0,600,600))

# par3d(userMatrix = view)
plot3d(x = river.ext[,2], y = river.ext[,3], z = river.ext[,1],
       xlab = "", ylab = "", zlab = "",
       axes = FALSE,cex=1.9)
arc3d(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0), 
                                  radius = 1, lwd = 2, col = "black")
arc3d(c(1, 0, 0), c(0, 0, 1), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(0, 0, 1), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
aspect3d(1, 1, 1)
axes3d(edges = c("x-+", "y-+", "z+-"),
       ntick = 6,                       # Attempt 6 tick marks on each side
       cex = 1)                       # Smaller font
# Add axis labels. 'line' specifies how far to set the label from the axis.
mtext3d("Station 28",       edge = "x-+", line = 5)
mtext3d("Station 29",       edge = "y-+", line = 2,pos=c(0,0,1))
mtext3d("Station 27",          edge = "z+-", line = 5)

pp <- dget("perspective.R")
par3d(pp)

##########################################
# K-means clustering on river extremes
##########################################

set.seed(800)
library(skmeans)
k <- 3
fit <- skmeans(river.ext,k,method="pclust",control = list(nruns = 1000, maxchains=100))
fit$cluster

open3d(windowRect = 0 + c( 0, 0, 600,600 ) )
plot3d(x = river.ext[fit$cluster==1,2], 
       y = river.ext[fit$cluster==1,3], 
       z = river.ext[fit$cluster==1,1],
       xlab = "", ylab = "", zlab = "",
       axes = FALSE,col='red',cex=1.9,pch=0)
plot3d(x = river.ext[fit$cluster==2,2], 
       y = river.ext[fit$cluster==2,3], 
       z = river.ext[fit$cluster==2,1],cex=1.9,pch=1,add=T)
plot3d(x = river.ext[fit$cluster==3,2], 
       y = river.ext[fit$cluster==3,3], 
       z = river.ext[fit$cluster==3,1],col='blue',cex=1.9,pch=2,add=T)
arc3d(c(1, 0, 0), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(1, 0, 0), c(0, 0, 1), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
arc3d(c(0, 0, 1), c(0, 1, 0), c(0, 0, 0), radius = 1, lwd = 2, col = "black")
aspect3d(1, 1, 1)
legend3d("topright", c("cluster 1", "cluster 2",'cluster 3'), 
         pch=16,col=c('black','red','blue'))
axes3d(edges = c("x-+", "y-+", "z+-"),
       ntick = 6,                       # Attempt 6 tick marks on each side
       cex = 1)                       # Smaller font
# Add axis labels. 'line' specifies how far to set the label from the axis.
mtext3d("Station 28",       edge = "x-+", line = 5)
mtext3d("Station 29",       edge = "y-+", line = 2,pos=c(0,0,1))
mtext3d("Station 27",          edge = "z+-", line = 5)

pp <- dget("perspective.R")
par3d(pp)