require(mclust)
require(MASS)
require(NbClust)
require(corrplot)
require(gplots)

# Local ~home directory
#home <- path.expand("~")
home <- "C:\\1-hongxiu\\codes\\github"

# Load data
setwd(paste0(home,"\\bootstrap\\"))
source("functions.R")
source("BASC_func.R")

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

asd_msn_corrected <- as.matrix(read.csv('meanMS_regional_corrected_asd.csv', header = FALSE))
Npsy.matrix <- asd_msn_corrected

sub_lists <- read.csv('Used_asd_control.csv')
sub_lists <- sub_lists[1:138,c("FILE_ID")]

D <- dim(asd_msn_corrected)

sample_idx <- sample(1:D[1], D[1], replace = TRUE)
M <- asd_msn_corrected[sample_idx,]

nc <- sample(2:30, 1) # create an integer randomly as the number of clusters
# clus <- MyNbClust_hx(M, nc)
partition <- kmeans(asd_msn_corrected, nc)#$cluster


organized_clusCut <- integer(D[1]) # create a vector whose values are all 0.
for (i in 1:length(sample_idx))
{
  organized_clusCut[sample_idx[i]] <- partition$cluster[i] # how to solve the repeated subjects are divided into different cluster?
  
}
clusCut <- organized_clusCut

# Obtains the row's IDs (removes X*.?) 
names(clusCut) <- sub_lists

# -----------------------------------------------------------------------
# PARAMETERS
Matrix <- asd_msn_corrected
boots <- 10

D <- dim(Matrix)   		# Dimensions of the matrix
Boots <- boots              	# N bootstraps
l <- N <- 0               	# l=loop counter & N=selected cluster partitions
k <- c()                 	# vector of best number of clusters selected by NbCluster

# Empty matrix (S0)
S0 <- (matrix(data=0,nrow = D[1],ncol = D[1])) # makes and empty matrix
colnames(S0) <- rownames(S0) <- sub_lists#rows.ID

# Total Sum of Stability matrices
S.sum <- S0


Sij.B <- mtx.jointProb(clusCut,S0) # Similarity boot matrix
S.bin <- ceiling(Sij.B/max(Sij.B)) # Binary ocurrence = Stability Mtx per boot
S.sum <- S.sum+S.bin # Sum of all Sij matrices



X <- clusCut
joint <- S0
J<-names(X)
for (i in names(X)) { 
  for (j in J) { 
    # Match the ID in the joint matrix
    x<-match(c(i,j),names(joint[1,]))
    if ( X[i]==X[j] ) { ij<-1 } else {ij<-0}
    val<-joint[x[1], x[2]]+ij
    joint[x[1], x[2]] <- joint[x[2], x[1]] <- val
    }
    
  
  J<-J[-match(i,J)]}



