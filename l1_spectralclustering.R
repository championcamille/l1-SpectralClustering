# l1_spectralclustering(): main code for running l1-spectral clustering

l1_spectralclustering <- function(M, k = NULL, index = NULL, pen, stab=TRUE){
  # M: the matrix we are working on (adjacency matrix, e.g. from CreateDataSet())
  # k: true number of clusters (not necessary needed). If not precised, k is chosen by spectral eigengap
  # index: representative elements of the clusters (not necessary needed). If not precised, index are chosen using the betweeness centrality score.
  # pen: penalty (to be chosen among "lasso" or "thresholdedLS")
  # stab: TRUE/FALSE, should the indices be stabilized?
  
  # outputs: 
  #   - comm: the community matrix
  #   - Structure: structure of the graph
  #   - clusters: number of clusters
  #   - Elements: the representative elements of the clusters
  
  # 1st step: finding the connected components
  Structure <- FindStructure(M)
  
  # 2nd step: finding the optimal number of clusters (only if k is not provided)
  clusters <- FindNbrClusters(M, Structure  = Structure, k = k)
  
  # 3rd step: finding the representative elements of the clusters
  Elements <- FindElement(M = M, Structure = Structure, clusters = clusters, index = index)
  
  # 4th step: running the l1-spectral clustering algorithm on each connected component (each of them are treated independtly)
  comm <- matrix(0,nrow=ncol(M),ncol=clusters$nbr_clusters_total)
  S <- cumsum(unlist(clusters$nbr_clusters)[-length(clusters$nbr_clusters)])
  for (i in (1:length(Structure$groups))){
    Mtmp <- M[Structure$groups[[i]],Structure$groups[[i]]]
    clusters_tmp <- clusters$nbr_clusters[[i]]
    indices_tmp <- Elements$indices[[i]][which(Elements$indices[[i]]%in%Structure$groups[[i]])]
    indices_tmp <- match(indices_tmp,Structure$groups[[i]])
    score_tmp <- Elements$score[[i]]
    names(score_tmp) <- paste0("Node",match(as.numeric(substring(names(Elements$score[[i]]),5)),Structure$groups[[i]]))
    Elements_tmp <- list(score = score_tmp,indices = indices_tmp)
    
    results <- l1spectral(M = Mtmp, k = clusters_tmp, elements = Elements_tmp, pen = pen, stab=stab)

    if (!is.null(ncol(results))){
      if (ncol(results)!=clusters$nbr_clusters[[i]]){
        clusters$nbr_clusters[[i]] <- ncol(results)
        S <- cumsum(unlist(clusters$nbr_clusters))
      }
    } else {
      if (clusters$nbr_clusters[[i]]!=1){
        clusters$nbr_clusters[[i]] <- 1
        S <- cumsum(unlist(clusters$nbr_clusters))
      } 
    }
    if (i==1){
      comm[Structure$groups[[i]],1:S[1]] <- results
    } else {
      comm[Structure$groups[[i]],(S[i-1]+1):S[i]] <- results
    }
  }
  
  if (length(which(colSums(comm)==0))>0){
    comm <- comm[,-which(colSums(comm)==0)]
  }
  comm[comm>0] <- 1
  return(list(comm=comm,Structure=Structure,clusters=clusters,Elements =Elements))
}

### function to find the connected components in a graph
FindStructure <- function(M){
  # M: adjacency matrix of the graph
  # Output: list of two elements
  #   - graph: igraph object with the graph
  #   - groups: list of connected components and corresponding nodes
  
  # define the graph
  graph <- graph_from_adjacency_matrix(M,mode="undirected")
  
  # find the connected components
  clu <- components(graph)
  groups <- groups(clu)
  if (length(groups)>1){
    groups <- groups[order(sapply(groups,length))]
  } else {
    groups <- groups
  }
  nbr_comp <- length(groups)
  names(groups) <- paste0("Component",c(1:nbr_comp))
  
  if (nbr_comp==1){
    print("The graph has only one connected component.")
  } else {
    print(paste0("The graph have ",nbr_comp," connected components. Each of them will be clustered using the l1-spectral clustering."))
  }
  return(list(graph=graph,groups=groups))
}

### function to estimate the optimal number of clusters (only if k=NULL)
FindNbrClusters <- function(M, Structure, k = NULL){
  # M: the adjacency matrix of the graph
  # Structure: already existing connected components 
  # k: number of clusters 
  # Outputs:
  #   - nbr_clusters: optimal number of clusters by component
  #   - nbr_clusters_total: optimal total number of clusters
  
  nbr_group <- length(Structure$groups) # number of connected components
  
  # various possible cases
  if (is.null(k)){
    # k is not provided
    
    eigenvalues <- lapply(X = c(Structure$groups,list(all=c(1:ncol(M)))),FUN = Eigen_list)
    
    par(mfrow=c(1,1))
    plot(eigenvalues$all,main="Eigenvalues of the Laplacian matrix",ylab="Eigenvalues",xlab="",type="b")
    
    if (nbr_group>1){
      for (i in (1:nbr_group)){
        points(eigenvalues[[i]],col=rainbow(nbr_group)[i],type="b")
      }
      legend("bottomright",legend = c("All nodes",paste0("Connected component ",c(1:nbr_group))),col=c("black",rainbow(nbr_group)[1:nbr_group]),lty=1)
    }
    
    gaps <- lapply(X = eigenvalues,FUN = Gap)
    par(xpd=FALSE) 
    for (i in (1:length(gaps))){
      color <- c(rainbow(nbr_group),"black")
      abline(v=gaps[[i]],col=color[i],ylim=c(0,10),lty=3)
    }
    
    if (length(gaps)==2){
      nbr_clusters_total <- gaps$all
      print(paste0("The optimal number of clusters is ",gaps$all,"."))
    } else {
      nbr_clusters_total <- sum(unlist(gaps)[1:(length(gaps)-1)])
      print(paste0("The optimal number of clusters is ",nbr_clusters_total,"."))
      print(paste0(c("Here,",nbr_group,"connected components were detected. Each of them should be clustered into",unlist(gaps[1:nbr_group]),"clusters."),collapse=" "))
    }
    
  } else if (nbr_group==1 && !is.null(k)){
    # easy case: one component and the number of clusters k is provided
    
    nbr_clusters_total <- k
    gaps <- list(Component1=k,all=k)
    print(paste0("The provided number of clusters is ",nbr_clusters_total,"."))
    
  } else {
    # last case: more than one component and the number of clusters k is provided

    if (k<nbr_group){
      
      print(paste0("The provided number of clusters ",k," is smaller than the number of components. This is not possible, we thus adjust the number of clusters to the number of components."))
      k <- nbr_group
      nbr_clusters_total <- k
      gaps <- rep(list(1),nbr_group)
      names(gaps) <- paste0("Component",1:nbr_group)
      gaps <- c(gaps, list(all=k))
      print(paste0("The total number of clusters is ",k,"."))
      
    } else if (k==nbr_group){
      
      k <- nbr_group
      nbr_clusters_total <- k
      gaps <- rep(list(1),nbr_group)
      names(gaps) <- paste0("Component",1:nbr_group)
      gaps <- c(gaps, list(all=k))
      print(paste0("The provided number of clusters is ",k,"."))
      
    } else {
      
      eigenvalues <- lapply(X = c(Structure$groups,list(all=c(1:ncol(M)))),FUN = Eigen_list)
      par(mfrow=c(1,1))
      plot(eigenvalues$all,main="Eigenvalues of the Laplacian matrix",ylab="Eigenvalues",xlab="",type="b")
      for (i in (1:nbr_group)){
        points(eigenvalues[[i]],col=rainbow(nbr_group)[i],type="b")
      }
      legend("bottomright",legend = c("All nodes",paste0("Connected component ",c(1:nbr_group))),col=c("black",rainbow(nbr_group)[1:nbr_group]),lty=1)
    
      gaps <- lapply(X = eigenvalues,FUN = Gap)
      par(xpd=FALSE) 
      for (i in (1:length(gaps))){
        color <- c(rainbow(nbr_group),"black")
        abline(v=gaps[[i]],col=color[i],ylim=c(0,10),lty=3)
      }
    
      if (sum(unlist(gaps)[1:(length(gaps)-1)]) != k){
        # there is a problem
        message(paste0("The optimal number of clusters ",sum(unlist(gaps)[1:(length(gaps)-1)])," does not coincide to the provided number of clusters."))
        gaps[[(length(gaps)-1)]] <- k-sum(unlist(gaps)[1:(length(gaps)-2)])
      }
      nbr_clusters_total <- sum(unlist(gaps)[1:(length(gaps)-1)])
      print(paste0("The provided number of clusters is ",k,"."))
      print(paste0(c("Here,",nbr_group,"connected components were detected. Each of them should be clustered into",unlist(gaps[1:nbr_group]),"clusters."),collapse=" "))
  
    }
  
  }

  return(list(nbr_clusters=gaps,nbr_clusters_total=nbr_clusters_total))
}

### function to find the representative elements 
FindElement <- function(M, Structure, clusters, index = NULL){
  # M: the adjacency matrix of the graph
  # Structure: Structure: already existing connected components 
  # clusters: output of the function FindNbrClusters()
  # index: representative elements of the clusters (not needed)
  # Outputs:
  #   - score: the edge betweenness score of all nodes
  #   - Nodes: list of representative elements
  
  n <- ncol(M)
  
  betweenness <- lapply(Structure$groups,between)
  
  if (is.null(index)){
    # case 1: no representative elements
    
    Nodes <- c()
    for (i in (1:length(betweenness))){ 
      b <- betweenness[[i]][1:unlist(clusters$nbr_clusters)[i]]
      b <- as.numeric(substring(names(rev(b)), 5))
      Nodes <- c(Nodes,list(b))
    }
    names(Nodes) <- names(Structure$groups)
  
  } else {
    # case 2: provided representative elements
    if (length(index)<clusters$nbr_clusters_total){
      print("The number of representative elements does not coincide with the number of clusters. Please add some elements to re-run the code. By default, new elements will be chosen.")
      
      Nodes <- c()
      for (i in (1:length(betweenness))){ 
        b <- betweenness[[i]][1:unlist(clusters$nbr_clusters)[i]]
        b <- as.numeric(substring(names(rev(b)), 5))
        Nodes <- c(Nodes,list(b))
      }
      names(Nodes) <- names(Structure$groups)
    } else{
      test <- lapply(Structure$groups,function(x){
             length(which(index%in%x))
      })
      test2 <- unlist(test) - unlist(clusters$nbr_clusters[-length(clusters$nbr_clusters)])
      
      if (length(which(test2<0))>0){
        print("At least one of the connected components does not have enough representative element. New elements will be chosen.")
        Nodes <- c()
        for (i in (1:length(betweenness))){ 
          b <- betweenness[[i]][1:unlist(clusters$nbr_clusters)[i]]
          b <- as.numeric(substring(names(rev(b)), 5))
          Nodes <- c(Nodes,list(b))
        }
        names(Nodes) <- names(Structure$groups)
      } else {
        Nodes <- c()
        for (i in (1:length(betweenness))){
          index_gr <- index[which(index %in% Structure$groups[[i]])]
          b <- betweenness[[i]][paste0("Node",index_gr)]
          if (length(b)>1){
            b <- sort(b,decreasing = TRUE)
          }
          if (test2[i] > 0){
            b <- b[1:clusters$nbr_clusters[[i]]]
          }
          b <- as.numeric(substring(names(rev(b)), 5))
          Nodes <- c(Nodes,list(b))
        }
        names(Nodes) <- names(Structure$groups)
      }
    }
  }
  return(list(score=betweenness,indices=Nodes))
}

### function to run the l1spectral on one component
l1spectral <- function(M, k, elements, pen, stab = TRUE){
  # M: the adjacency matrix of the graph
  # k: the number of clusters to form (output of the function FindClusters())
  # elements: representative elements of the clusters (output of the function FindElement())
  # pen: the penalty (to be chosen among lasso and threshold)
  # stab: should the indices stabilized? TIME CONSUMING
  # Outputs:
  #   comm: matrix of the components
  
  if (length(M)==1){
    # only one node in the community
    comm <- 1
  } else {
    # code for running the l1-spectral algorithm for one component
    indices <- elements$indices
    
    # 1st step: svd on M
    n <- ncol(M)
    svd <- eigen(M) 
    eigenvalues <- sort(svd$values,index.return=TRUE)
    eigenvectors <- svd$vectors[,eigenvalues$ix]
    
    # 2nd step: loop on the number of clusters  
    algo <- "stop"
    comm <- c()
    DoubleNodes <- c()
    
    while (algo == "stop"){
      if (length(DoubleNodes)>0){
        # find other indices
        I <- elements$score[-which(names(elements$score)%in%DoubleNodes)]
        if (length(I)==0){
          print("One cluster disappears.")
          DoubleNodes <- c()
          algo <- "continue"
          break
        } else if (length(I)<k){
          print("One cluster disappears.")
          DoubleNodes <- c()
          comm <- c()
          k <- k-1
          indices <- elements$indices[2:(k+1)]
        } else {
          I <- names(sort(I[1:k]))
          indices <- as.numeric(substring(I, 5))
        }
      }
    
      if (k>1){
        eigenvectors_tmp <- eigenvectors
        comm <- c()
        for (i in (1:k)){
          # 3rd step: check the indices (only if i>1)
          if (stab==TRUE){
            if (i>1){
              if (length(which(v[indices[-(i-1)]]>0))>0){
                print("Find other community indices.")
                doubleNodes <- paste0("Node",indices[i-1])
                DoubleNodes <- c(DoubleNodes,doubleNodes)
                algo <- "stop"
                break
              } else {
                algo <- "continue"
              }
            } 
          }
        
          # 4th step: Gram-Schmidt (only if i>1)
          if (i>1){
            eigenvectors_tmp <- eigenvectors_tmp[,-(n-k+i-1)]
            eigenvectors_tmp <- cbind(v,eigenvectors_tmp)
          
           eigenvectors_tmp <- GramSchmidt(eigenvectors_tmp)
           eigenvectors_tmp <- cbind(eigenvectors_tmp[,2:(n-k+i-1)],eigenvectors_tmp[,1],eigenvectors_tmp[,(n-k+i):n])
          }
        
          # 5th step: solve the lasso
          U <- t(eigenvectors_tmp[,1:(n-k+i-1)])
          v <- Lasso(U, n, indices, iteration = i, pen=pen, k)
          print(paste0("Iteration ",i," done."))
        
          # 6th step: save the community index
          comm <- cbind(comm,v)
        
          I <- which(rowSums(comm>0)>1)
          if (length(I)>0){
            for (j in (1:length(I))){
              C <- comm[I[j],] 
              C[C<max(C)] <- 0
              comm[I[j],] <- C
            }
          } 
          if (stab==TRUE){
            if (i==k){
              # check the indices for the last time
              if (length(which(v[indices[-i]]>0))>0){
                print("Find other community indices.")
                doubleNodes <- paste0("Node",indices[-i][which(v[indices[-i]]>0)])
                DoubleNodes <- c(DoubleNodes,doubleNodes)
                algo <- "stop"
                break
              } else {
                algo <- "continue"
              }
            }
          } else {
            algo <- "continue"
          }
        }
      } else {
        comm <- rep(1,n)
        algo <- "continue"
      }
    }
  }
  return(comm)
}

# useful functions
Eigen_list <- function(group){
  # spectral decomposition on subgraphs of a graph
  if (length(group)>1){
    D <- diag(rowSums(M[group,group])) # degree matrix
    L <- D-M[group,group]
    svd <- eigen(L) 
    eigenvalues <- sort(svd$values)
  } else {
    eigenvalues <- NA
  }
  return(eigenvalues)
}

Gap <- function(eigenvalue){
  # compute the spectral gap
  if (length(eigenvalue)==1){
    nbr_cluster <- 1
  } else {
    gap <- c()
    for (i in (2:length(eigenvalue))){
      gap <- c(gap,eigenvalue[i]-eigenvalue[i-1])
    }
    gap_ecart <- c(gap,0)-c(0,gap)
    nbr_cluster <- which(gap_ecart[2:length(gap_ecart)]>0.20)[1]+1
    if (is.na(nbr_cluster)){
      nbr_cluster=1
    }
  }
  return(nbr_cluster)
}

between <- function(comm){
  # betweeness centrality score on subgraphs of a graph
  graph_small <- graph_from_adjacency_matrix(M[comm,comm],mode="undirected")
  b <- betweenness(graph_small,directed = FALSE,normalized=TRUE)
  names(b) <- paste0("Node",comm)
  I <- order(b,decreasing = TRUE)
  b <- b[I]
  return(b)
}

Lasso <- function(U, n, indices, iteration,pen,k){
  # solve the l1-min problem using Lasso or thresholded least-squares penalization
  w <- U[,indices[iteration]]
  W <- matrix(U[,-indices[iteration]],ncol=(ncol(U)-1),nrow=nrow(U))
  
  if (sum(w)==0 || length(w)==1){
    print("There is only one node in this cluster.")
    # w is constant - no more nodes in the cluster
    v <- rep(0,n)
    v[indices] <- 1
  } else if (length(which(w!=0))==1){
    sol <- glmnet(W,-w,lambda=0,lower.limits=0)
    sol <- sol$beta
    solution <- as.matrix(sol)
    
    solution[solution<0.5] <- 0
    
    if(indices[iteration]==1){
      v <- c( 1 ,solution[indices[iteration]:length(solution)])
    } else{
      v <- c(solution[1:(indices[iteration]-1)], 1 ,solution[indices[iteration]:length(solution)])
    }
    
  } else {
    if (pen == "lasso"){
      lassosol <- cv.glmnet(W,-w)
      
      cvtop<- min(lassosol$cvm)+5*(max(lassosol$cvm)-min(lassosol$cvm))/100
      plot(lassosol)
      abline(h=cvtop)
      
      error <- min(lassosol$cvm[lassosol$cvm>cvtop])
      lambda_opt <- lassosol$lambda[which(lassosol$cvm==error)]
      
      sol <- glmnet(W,-w,lambda=lambda_opt)
      sol <- sol$beta
      solution <- as.matrix(sol)
      
      I <- which(solution!=0)
      sol2 <- glmnet(W[,I],-w,lambda=0)
      sol2 <- sol2$beta
      sol2 <- as.matrix(sol2)
      solution[I] <- sol2
    } else {
      # thresholded least-squares penalty
      # cross validation test
      groups <- cvFolds(nrow(W),K=5)
      T <- seq(from=0, to=0.5,0.001)
      erreur_cv <- c()
      for (g in (1:5)){
        Erreur <- c()
        sol <- glmnet(W[which(groups$which!=g),],-w[which(groups$which!=g)],lambda=0)
        sol <- sol$beta
        sol <- as.matrix(sol)
        
        for (t in (1:length(T))){
          sol2 <- sol
          sol2[which(sol<T[t])] <- 0
          if (length(which(sol>0))>1){
            solution <- glmnet(W[which(groups$which!=g),which(sol>0)],-w[which(groups$which!=g)],lambda=0)
            solution <- solution$beta
          } else {
            solution <- lm(-w[which(groups$which!=g)]~W[which(groups$which!=g),which(sol>0)]-1)
            solution <- solution$coefficients
          }
          sol2[which(sol2>0)] <- solution
          
          erreur <- sum((-w[which(groups$which==g)]-W[which(groups$which==g),]%*%sol2)^2)
          Erreur <- c(Erreur,erreur)  
        }
        
        erreur_cv <- cbind(erreur_cv,Erreur)
      }
      Erreur_cv <- rowMeans(erreur_cv)
      
      plot(T,Erreur_cv,type="l")
      lambda_opt <- T[which.min(Erreur_cv)]
      sol2 <- glmnet(W,-w,lambda=0)
      sol2 <- sol2$beta
      sol2 <- as.matrix(sol2)
      sol2[which(sol2<lambda_opt)] <- 0
      solution <- sol2
    }
    if(indices[iteration]==1){
      v <- c( 1 ,solution[indices[iteration]:length(solution)])
    } else if (indices[iteration]==n){
      v <- c(solution[1:(indices[iteration]-1)], 1)
    } else {
      v <- c(solution[1:(indices[iteration]-1)], 1 ,solution[indices[iteration]:length(solution)])
    }
  }
  v[v<0] <- 0
  return(v)
}

GramSchmidt <- function(U){
  # Gram-Schmidt function
  B <- grahm_schimdtCpp(U)
  B <- B$Q
  return(B)
}

