##################################
# This script contains the functions for applying spherical K-means,  spherical k-pc, and spectral clustering for extremes.
# All algorithms are implemented on angular components estimated from the data.  That is, data = angular components.
##################################

##################################
# Spherical K-Means function

library(skmeans)
clusterMeans <- function(data,k=2,nruns=500){
  fit <- skmeans(data,k,method="pclust",control=list(nruns = nruns))
  return(list(centers=fit$prototypes,assignments=fit$cluster))
} 


##################################
# Spectral clustering function

clusterSpectral <- function(data,k=2,k.neighbour=10){
  dist_matrix <- as.matrix(dist(data,diag=T,upper=T))
  similarity_matrix <- exp(-dist_matrix^2/(2 * 1^2))
  n.ext <- nrow(data)
  W <- matrix(0,n.ext,n.ext)
  for (i in 1:n.ext){
    sorted.index <- order(similarity_matrix[i,],decreasing=T)
    keep.index <- sorted.index[1:k.neighbour+1]
    W[keep.index,i] <- similarity_matrix[keep.index,i]
    W[i,keep.index] <- similarity_matrix[i,keep.index]
  }
  # W <- similarity_matrix
  D12 <- diag(1/sqrt(rowSums(W)))
  L <- diag(1,nrow=n.ext) - D12 %*% W %*% D12
  eigenvectors <- eigen(L)$vectors
  U <- eigenvectors[, -(1:(n.ext-k))]
  for (i in 1:n.ext){
    U[i,] <- U[i,]/sqrt(sum(U[i,]^2))
  }
  cluster_assignments <- kmeans(U, centers = k)$cluster
  return(cluster_assignments)
} 

# Calculate cluster centers from cluster assignments

clusterCenters <- function(data,cluster_assignments){
  n.centers <- length(unique(cluster_assignments))
  centers <- matrix(NA,n.centers,4)
  for (m in 1:n.centers){
    if (is.vector(data[cluster_assignments==m,]) == TRUE){
      centers[m,] <- data[cluster_assignments==m,]
    }
    else {
      centers[m,] <- colMeans(data[cluster_assignments==m,])
    }
  }
  return(centers)
}


#####################################################################
# Spherical K-Pc clustering
# The following code is taken from the supplementary material of Fomichov and Ivanov (2021)

library(ramify)
#######################################################################
################## spherical k-principal component clustering##########
#The main routine:
#clusterPC=function(data,k=2,nrep=100,tol=10^(-5),startFromMeans=FALSE)
#
#data is a matrix of observations, where each row is non-negative with unit Euclidean norm.
#k is the number of clusters
#nrep is a number of random restarts
#tol is used to stop searching for a local optimum for each restart
#startFromMeans=TRUE adds one restart using k-means centroids

#######################
#single iteration
#centroids is a k*d matrix with current proposals
clusterPC_iter=function(data, centroids){
  k=length(centroids[,1])
  n=length(data[,1])
  d=length(data[1,])
  M=data%*%t(centroids)
  #find current value
  v=mean(apply(M,1,max))
  gr=argmax(M,rows=T)
  for (i in 1:k){
    seldata=data[gr==i,]
    if (length(seldata)==d) #interpretation problem when just one vector
      seldata=t(seldata)
    Sig=t(seldata)%*%seldata/n
    res=eigen(Sig)
    centroids[i,]=abs(res$vectors[,1]) #use the first eigenvector, the entries fo which are necessarily positive
  }
  list(centroids,v)
}
#pick randomly the initial centers
clusterPCOnce=function(data,k,tol,startFromMeans=FALSE){
  val=0
  n=length(data[,1])
  if (startFromMeans)
    centroids=clusterMeans(data,k)  
  else{
    centroids=data[sample(1:n,k),]
    if (k==1)   #make sure it is a matrix
      centroids=t(as.matrix(centroids))
  }
  niter=0
  repeat{
    niter=niter+1
    res=clusterPC_iter(data, centroids)
    centroids=res[[1]]
    diff=res[[2]]-val
    val=res[[2]]
    if(diff<tol)
      break
  }
  #print(niter)
  list(centroids,val)
}
#iterate nrep times and pick the best
clusterPC=function(data,k=2,tol=10^(-5),nrep=100,startFromMeans=FALSE){
  maxval=0
  for( i in 1:nrep){
    res=clusterPCOnce(data,k,tol,startFromMeans && (i==1))
    if (res[[2]]>maxval){
      maxval=res[[2]]
      centroids=res[[1]]
    }
  }
  centroids
}

####################Useful routines#################
#assign groups
getClusterIndex=function(data,centroids){
  M=data%*%t(centroids)
  gr=argmax(M,rows=T)
  gr
}
#compute dissimilarity cost
getCost=function(data,centroids,cosine=TRUE){
  M=data%*%t(centroids)
  products=apply(M,1,max)
  if(cosine)
    return(1-mean(products))
  else
    return(1-mean(products^2))
}

