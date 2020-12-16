# load the libraries
library(igraph)
library(cvTools)
library(NMI)

# create a dataset
Data <- CreateDataSet(k=3,n=20,p=list(p_inside=0.1,p_outside=0.1))
M <- Data$A_hat

# load a toy dataset
load("ToyData.Rdata")
M <- Data$A_hat

# find the structure of the graph
Structure <- FindStructure(M)

# find the optimal number of clusters
clusters <- FindNbrClusters(M, Structure  = Structure)

# find the indices of the components
Elements <- FindElement(M = M, Structure = Structure, clusters = clusters)
Elements

# run the l1-spectral clustering
results <- l1_spectralclustering(M = M, pen = "thresholdedLS")
results$comm

# or unstabilized version
results2 <- l1_spectralclustering(M = M, pen = "thresholdedLS", stab=FALSE)
results2$comm

# with predefined clusters and indices
results3 <- l1_spectralclustering(M = M, pen = "thresholdedLS", k=4, index = c(13,2,7,15))
results3$comm

# compute the performances
Perfs <- ComputePerformances(Results=results$comm, A=Data$A)
Perfs
