####################################
# Example in Section 12.6.1
####################################

## Load algorithm functions
source('HelperFunctions.R')

## Simulate Data
set.seed(100)

d <- 4

b1 <- c(1,0,0,0)
b2 <- c(0,38/40,1/40,1/40)
b3 <- c(0,1/40,38/40,1/40)
b4 <- c(0,1/40,1/40,38/40)

# Function for detecting whether or not the estimated centers reside on the destinated subfaces
subfaces.detect <- function(centers){
  n.centers <- nrow(centers)
  is.subface1 <- c()
  is.subface2 <- c()
  for (m in 1:n.centers){
    center <- centers[m,]
    is.subface1[m] <- center[1]>0 & sum(center[-1])==0
    is.subface2[m] <- center[1]==0 & prod(center[-1])>0
  }
  return(sum(is.subface1)*sum(is.subface2)>0)
}


#############################

# First, q = 0.95

q <- 0.95

detection.rate1 <- c()
detection.rate2 <- c()
detection.rate3 <- c()

for (j in 1:100){
  n <- 10000
  Z1 <- (-log(runif(n)))^(-1)
  Z2 <- (-log(runif(n)))^(-1)
  Z3 <- (-log(runif(n)))^(-1)
  Z4 <- (-log(runif(n)))^(-1)
  X <- matrix(NA,n,d) 
  for (i in 1:n){
    X[i,] <- apply(rbind(b1*Z1[i],b2*Z2[i],b3*Z3[i],b4*Z4[i]),2,sum)
  }
  
  # Process data
  # Transform and truncate
  Paretotrans<-function(x) 1/(1-ecdf(x)(x)*length(x)/(length(x)+1)) 
  X.ext <- apply(X,2,Paretotrans)
  norm_vec <- function(x) sqrt(sum(x^2)) 
  norms<-apply(X.ext,1,norm_vec)
  X.ext <- X.ext[norms>quantile(norms,q),]
  X.ext<- X.ext/norms[norms>quantile(norms,q)]
  
  ###########################
  # k-means
  ###########################
  # set.seed(800)
  library(skmeans)
  k <- 2
  centers <- clusterMeans(X.ext,k)$centers
  centers[centers<0.15] <- 0
  hit <- subfaces.detect(centers)
  detection.rate1 <- c(detection.rate1,hit)
  
  ###########################
  # k-pc
  ###########################
  # source('KpcFns.R')
  centers <- clusterPC(X.ext,k)
  
  centers[centers<0.1] <- 0
  hit <- subfaces.detect(centers)
  detection.rate2 <- c(detection.rate2,hit)
  
  ###########################
  # spectral clustering
  ###########################
  
  # set.seed(800)
  cluster_assignments <- clusterSpectral(X.ext,k,20)
  
  centers <- clusterCenters(X.ext,cluster_assignments)
  
  centers[centers<0.1] <- 0
  hit <- subfaces.detect(centers)
  detection.rate3 <- c(detection.rate3,hit)
  # print(centers)
  print(j)
}

sum(detection.rate1)/100 
sum(detection.rate2)/100 
sum(detection.rate3)/100 

##################################
# Repeat for q=0.995

## Simulate Data

set.seed(100)

q <- 0.995

detection.rate1 <- c()
detection.rate2 <- c()
detection.rate3 <- c()

for (j in 1:100){
  
  n <- 10000
  Z1 <- (-log(runif(n)))^(-1)
  Z2 <- (-log(runif(n)))^(-1)
  Z3 <- (-log(runif(n)))^(-1)
  Z4 <- (-log(runif(n)))^(-1)
  X <- matrix(NA,n,d) 
  for (i in 1:n){
    X[i,] <- apply(rbind(b1*Z1[i],b2*Z2[i],b3*Z3[i],b4*Z4[i]),2,sum)
  }
  
  
  # Process data
  
  # Transform and truncate
  Paretotrans<-function(x) 1/(1-ecdf(x)(x)*length(x)/(length(x)+1)) 
  X.ext <- apply(X,2,Paretotrans)
  norm_vec <- function(x) sqrt(sum(x^2)) 
  norms<-apply(X.ext,1,norm_vec)
  
  q <- 0.995
  X.ext <- X.ext[norms>quantile(norms,q),]
  X.ext<- X.ext/norms[norms>quantile(norms,q)]
  
  ###########################
  # k-means
  ###########################
  # set.seed(800)
  library(skmeans)
  k <- 2
  centers <- clusterMeans(X.ext,k)$centers
  centers[centers<0.1] <- 0
  hit <- subfaces.detect(centers)
  detection.rate1 <- c(detection.rate1,hit)
  
  ###########################
  # k-pc
  ###########################
  # source('KpcFns.R')
  centers <- clusterPC(X.ext,k)
  centers[centers<0.1] <- 0
  hit <- subfaces.detect(centers)
  detection.rate2 <- c(detection.rate2,hit)
  
  ###########################
  # spectral clustering
  ###########################
  
  # set.seed(800)
  cluster_assignments <- clusterSpectral(X.ext,k,20)
  
  centers <- clusterCenters(X.ext,cluster_assignments)
  
  centers[centers<0.1] <- 0
  hit <- subfaces.detect(centers)
  detection.rate3 <- c(detection.rate3,hit)
  
  print(j)
}

sum(detection.rate1)/100 
sum(detection.rate2)/100 
sum(detection.rate3)/100 
