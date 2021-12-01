# # # # ###################################################################################################### # # # #
# 				Boostrap Analysis of Stable Clusters  
# 		Stratified Bootstrap Sampling of a Vector (Y), & Generation of a S.ij^{boot} Matrix
# # # # ###################################################################################################### # # # #
# Raúl Rodríguez Cruces
# raulrcruces@inb.unam.mx
# Created on May 2017
# Modified: January 2018

BASC <- function(Matrix=NULL, rows.ID=NULL, boots=NULL, type="ward2.D2"){
	# Matrix	Observations x Variables matrix of standarized values (ex. z-scored)
	# row.ID	Vector of strings with the ID of each observarion (row)
	# boots		Number of bootstraps to be perform	
	# type    "waard2.D2" or random-K stochastic
  #
	# All arguments are mandatory
	if (is.null(Matrix)) 
        stop("[ERROR]... The Matrix is NULL")
	if (is.null(rows.ID)) 
        stop("[ERROR]... The row's ID vector is is NULL")
	if (is.null(boots)) 
        stop("[ERROR]... Number of bootstraps is EMPTY")

	# Require the modified version of NbClust.R that allows the option: plotetc = FALSE
	# source(P("/Users/labuser/Desktop/hongxiu/2020_cognition_connectomics_TLE/code/BASC/functions.R"))	
  source("C:\\1-hongxiu\\codes\\github\\bootstrap\\functions.R")

	# -----------------------------------------------------------------------
	# FUNCTION - Creates stability matrix from a clustering partition vector
	mtx.jointProb <- function(clustPartition, joint) {
	  X<-clustPartition
	  J<-names(X)
	  for (i in names(X)) { 
	    for (j in J) { 
	      # Match the ID in the joint matrix
	      x<-match(c(i,j),names(joint[1,]))
	      if ( X[i]==X[j] ) { ij<-1 } else {ij<-0}
	      val<-joint[x[1], x[2]]+ij
	      joint[x[1], x[2]] <- joint[x[2], x[1]] <- val}
	    J<-J[-match(i,J)]}
	  return(joint)}

	# -----------------------------------------------------------------------
	# PARAMETERS
	D <- dim(Matrix)   		# Dimensions of the matrix
  	Boots <- boots              	# N bootstraps
  	l <- N <- 0               	# l=loop counter & N=selected cluster partitions
  	k <- c()                 	# vector of best number of clusters selected by NbCluster
	
	# Empty matrix (S0)
	S0 <- (matrix(data=0,nrow = D[1],ncol = D[1])) # makes and empty matrix
	colnames(S0) <- rownames(S0) <- rows.ID

	# Total Sum of Stability matrices
	S.sum <- S0

	# -----------------------------------------------------------------------
	# BOOTSTRAP  
	  while ( mean(diag(S.sum)) < Boots ) {l=l+1 #iteration counter

		# Stratified Bootstrap of Y, with replacement
		M <- t(sample(data.frame(t(Matrix)),D[1], replace=TRUE))
		clus <- 0   # Is it a good clustering??
		
		try(

		# Clustering algorithm & best partition based on NbClust
		    clus<-MyNbClust(M, distance = "euclidean", min.nc=2, max.nc=D[2]-1, 
		                    method = "ward.D2", index = "all",plotetc = FALSE))
		if (class(clus) != "numeric" ){ N=N+1
		    # Obtains the best partition
		    clusCut <- clus$Best.partition
		    k <- c(k,max(clusCut))

		    # Obtains the row's IDs (removes X*.?) 
		    s <- names(clusCut)
		    s <- sapply(strsplit(s, split='X', fixed=TRUE), function(x) (x[2]))
		    names(clusCut) <- s <- sapply(strsplit(s, split='.', fixed=TRUE), function(x) (x[1]))

		    # Stability matrix
		    Sij.B <- mtx.jointProb(clusCut,S0) # Similarity boot matrix
		    S.bin <- ceiling(Sij.B/max(Sij.B)) # Binary ocurrence = Stability Mtx per boot
		    S.sum <- S.sum+S.bin # Sum of all Sij matrices
		}

		print("*************************************************")
		print(paste("Loop",l,", Selected partitions",N,"... mean ij occurrence:",mean(diag(S.sum))))
		print("*************************************************")

		
	  }

	# -----------------------------------------------------------------------
	# RESULTS  
	# Vector of i occurrence is the diagonal of the Total Sum of stability matrices (Sij)
	return(list(Sij=S.sum, k=k, N=N))
}


BASC_hx <- function(Matrix=NULL, rows.ID=NULL, boots=NULL, type="ward2.D2"){
  # Matrix	Observations x Variables matrix of standarized values (ex. z-scored)
  # row.ID	Vector of strings with the ID of each observarion (row)
  # boots		Number of bootstraps to be perform	
  # type    "waard2.D2" or random-K stochastic
  #
  # All arguments are mandatory
  if (is.null(Matrix)) 
    stop("[ERROR]... The Matrix is NULL")
  # if (is.null(rows.ID)) 
  #   stop("[ERROR]... The row's ID vector is is NULL")
  if (is.null(boots)) 
    stop("[ERROR]... Number of bootstraps is EMPTY")
  
  # Require the modified version of NbClust.R that allows the option: plotetc = FALSE
  # source(P("/Users/labuser/Desktop/hongxiu/2020_cognition_connectomics_TLE/code/BASC/functions.R"))	
  source("C:\\1-hongxiu\\codes\\github\\bootstrap\\functions.R")
  
  # -----------------------------------------------------------------------
  # FUNCTION - Creates stability matrix from a clustering partition vector
  mtx.jointProb <- function(clustPartition, joint) {
    X<-clustPartition
    J<-names(X)
    for (i in names(X)) { 
      for (j in J) { 
        # Match the ID in the joint matrix
        x<-match(c(i,j),names(joint[1,]))
        if ( X[i]==X[j] ) { ij<-1 } else {ij<-0}
        val<-joint[x[1], x[2]]+ij
        joint[x[1], x[2]] <- joint[x[2], x[1]] <- val}
      J<-J[-match(i,J)]}
    return(joint)}
  
  # -----------------------------------------------------------------------
  # PARAMETERS
  D <- dim(Matrix)   		# Dimensions of the matrix
  Boots <- boots              	# N bootstraps
  l <- N <- 0               	# l=loop counter & N=selected cluster partitions
  k <- c()                 	# vector of best number of clusters selected by NbCluster
  
  # Empty matrix (S0)
  S0 <- (matrix(data=0,nrow = D[1],ncol = D[1])) # makes and empty matrix
  colnames(S0) <- rownames(S0) <- rows.ID
  
  # Total Sum of Stability matrices
  S.sum <- S0
  
  # -----------------------------------------------------------------------
  # BOOTSTRAP  
  while ( mean(diag(S.sum)) < Boots ) {l=l+1 #iteration counter
  
  # Stratified Bootstrap of Y, with replacement
  M <- t(sample(data.frame(t(Matrix)),D[1], replace=TRUE))
  clus <- 0   # Is it a good clustering??
  
  # try(
    
    # Clustering algorithm & best partition based on NbClust
    # clus<-MyNbClust(M, distance = "euclidean", min.nc=2, max.nc=D[2]-1, 
    #                 method = "ward.D2", index = "all",plotetc = FALSE))
  nc <- sample(2:137, 1) # create an integer randomly as the number of clusters
  clus <- MyNbClust_hx(M, nc)
  
  #############
  # clus is the direct output of k-means, needs to be changed
  #############
    

  if (class(clus) != "numeric" ){ N=N+1
  # Obtains the best partition
  clusCut <- clus#$Best.partition
  k <- c(k,max(clusCut))
  
  # Obtains the row's IDs (removes X*.?) 
  s <- names(clusCut)
  s <- sapply(strsplit(s, split='X', fixed=TRUE), function(x) (x[2]))
  names(clusCut) <- s <- sapply(strsplit(s, split='.', fixed=TRUE), function(x) (x[1]))
  
  # Stability matrix
  Sij.B <- mtx.jointProb(clusCut,S0) # Similarity boot matrix
  S.bin <- ceiling(Sij.B/max(Sij.B)) # Binary ocurrence = Stability Mtx per boot
  S.sum <- S.sum+S.bin # Sum of all Sij matrices
  }
  
  print("*************************************************")
  print(paste("Loop",l,", Selected partitions",N,"... mean ij occurrence:",mean(diag(S.sum))))
  print("*************************************************")
  
  
  }
  
  # -----------------------------------------------------------------------
  # RESULTS  
  # Vector of i occurrence is the diagonal of the Total Sum of stability matrices (Sij)
  return(list(Sij=S.sum, k=k, N=N))
}







