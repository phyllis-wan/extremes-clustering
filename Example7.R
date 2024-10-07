###########################
# Example in Section 12.7
###########################

source('HelperFunctions.R')

# Load data
library(graphicalExtremes)
data <- danube$data_clustered
n <- nrow(data)
dat <- data
d <- ncol(data)

# Transform and truncate
Frechettrans<-function(x) 1/(1-ecdf(x)(x)*length(x)/(length(x)+1)) 
river.ext <- apply(dat,2,Frechettrans)
norm_vec <- function(x) sqrt(sum(x^2)) 
norms<-apply(river.ext,1,norm_vec)

q <- 0.9
river.ext<-river.ext[norms>quantile(norms,q),]
norms<-apply(river.ext,1,norm_vec)
river.ext<-river.ext/norms 


set.seed(200)
## Elbow plot with k-means ##
value<-rep(0,30)
for (k in 1:30){
  value[k]<-skmeans(river.ext,k,method="pclust",control = list(nruns = 500))$value
  print(k)
}
par(mfrow=c(1,1))
pdf('Elbow.pdf',10,6)
plot(1:30,value[1:30],ylab='objective function', xlab="K")
dev.off()


###########################
# k-means
###########################

set.seed(300)

library(skmeans)
k <- 6
fit <- clusterMeans(river.ext,k)
fit$assignment

# Re-order the clusters.  Here I simply manually re-labeled them.  Not the smartest move.
cluster.id <- rep(NA,length(fit$cluster))
cluster.id[fit$assignment==1] <- 2
cluster.id[fit$assignment==2] <- 6
cluster.id[fit$assignment==3] <- 1
cluster.id[fit$assignment==4] <- 5
cluster.id[fit$assignment==5] <- 4
cluster.id[fit$assignment==6] <- 3

# Plot the observations in heatmap
library('plot.matrix')
test.matrix <- river.ext[order(cluster.id),]
cluster <- sort(cluster.id)
row.names(test.matrix) <- as.character(sort(cluster.id))
colnames(test.matrix) <- as.character(1:31)
# pdf('1Kmeans.pdf',12,10)
plot(test.matrix,main='Spherical K-means clustering')
# dev.off()

test.matrix1 <- c()
rowname.vec <- c()
for (j in 1:6){
  test.matrix1 <- rbind(test.matrix1,river.ext[cluster.id==j,])
  test.matrix1 <- rbind(test.matrix1,rep(-2,d))
  rowname.vec <- c(rowname.vec,rep(as.character(j),sum(cluster.id==j)),'')
}
row.names(test.matrix1) <- rowname.vec 
colnames(test.matrix1) <- as.character(1:31)


pdf('1Kmeans.pdf',12,10)
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(test.matrix1,main='Spherical K-means clustering',breaks=range(river.ext),ylab='Clusters',xlab='Components')
dev.off()



 ###########################
# k-pc
###########################


set.seed(400)

# Use the KPC algorithm from Fomichov and Ivanov
fitPC <- clusterPC(river.ext,k)

# Extract cluster allocations
n.ext <- nrow(river.ext)
cluster.id.PC <- rep(NA,n.ext)
for (i in 1:n.ext){
  cluster.id.PC[i] <- which.max(fitPC%*%river.ext[i,])
}

# Re-order
cluster.id <- rep(NA,n.ext)
cluster.id[cluster.id.PC==1] <- 6
cluster.id[cluster.id.PC==2] <- 3
cluster.id[cluster.id.PC==3] <- 2
cluster.id[cluster.id.PC==4] <- 5
cluster.id[cluster.id.PC==5] <- 4
cluster.id[cluster.id.PC==6] <- 1

# Plot
# After each cluster, an empty row to distinguish the clusters
library('plot.matrix')
test.matrix <- river.ext[order(cluster.id),]
cluster <- sort(cluster.id)
row.names(test.matrix) <- as.character(sort(cluster.id))
colnames(test.matrix) <- as.character(1:31)
plot(test.matrix,main='Spherical K-pc clustering',breaks=range(river.ext),ylab='Clusters',xlab='Components')

test.matrix1 <- c()
rowname.vec <- c()
for (j in 1:6){
  test.matrix1 <- rbind(test.matrix1,river.ext[cluster.id==j,])
  test.matrix1 <- rbind(test.matrix1,rep(-2,d))
  rowname.vec <- c(rowname.vec,rep(as.character(j),sum(cluster.id==j)),'')
}
row.names(test.matrix1) <- rowname.vec 
colnames(test.matrix1) <- as.character(1:31)


pdf('1Kpc.pdf',12,10)
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(test.matrix1,main='Spherical K-means and K-pc clustering',breaks=range(river.ext),ylab='Clusters',xlab='Components')
dev.off()

###########################
# spectral clustering
###########################

set.seed(500)

cluster_assignments <- clusterSpectral(river.ext,k,10)

cluster.id <- rep(NA,n.ext)
cluster.id[cluster_assignments ==1] <- 1
cluster.id[cluster_assignments ==2] <- 6
cluster.id[cluster_assignments ==3] <- 2
cluster.id[cluster_assignments ==4] <- 4
cluster.id[cluster_assignments ==5] <- 3
cluster.id[cluster_assignments ==6] <- 5
library('plot.matrix')
test.matrix <- river.ext[order(cluster.id),]
cluster <- sort(cluster.id)
row.names(test.matrix) <- as.character(sort(cluster.id))
colnames(test.matrix) <- as.character(1:31)

plot(test.matrix,main='Spectral clustering',breaks=range(river.ext),ylab='Clusters',xlab='Components')


test.matrix1 <- c()
rowname.vec <- c()
for (j in 1:6){
  test.matrix1 <- rbind(test.matrix1,river.ext[cluster.id==j,])
  test.matrix1 <- rbind(test.matrix1,rep(-2,d))
  rowname.vec <- c(rowname.vec,rep(as.character(j),sum(cluster.id==j)),'')
}
row.names(test.matrix1) <- rowname.vec 

colnames(test.matrix1) <- as.character(1:31)


pdf('1Spectral.pdf',12,10)
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(test.matrix1,main='Spectral clustering',breaks=range(river.ext),ylab='Clusters',xlab='Components')
dev.off()
