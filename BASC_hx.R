#################################################################
#   Clustering with best number of clusters and stable analysis
#################################################################
# https://www.rdocumentation.org/packages/NbClust/versions/3.0/topics/NbClust
# https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html
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

asd_msn_corrected <- as.matrix(read.csv('meanMS_regional_corrected_asd.csv', header = FALSE))
Npsy.matrix <- asd_msn_corrected

sub_lists <- read.csv('Used_asd_control.csv')
sub_lists <- sub_lists[1:138,c("FILE_ID")]

# ----------------------------------------------------------------------------------- #
#       BASC - Boostrap Analysis of Stable Clusters
# ----------------------------------------------------------------------------------- #

k10 <- BASC_hx(Matrix = Npsy.matrix, rows.ID = sub_lists, boots = 10)
# 
# # Saves the result of the 10,000 bootstrap to a .RData file
# #save(k10,file = paste0(home,"/git_here/micasoft/sandbox/raul/BASC/10k_34subjs.RData"))
# 
# #Loads the data
# # load(paste0(home,"/git_here/micasoft/sandbox/raul/BASC/10k_34subjs_new.RData"))
# 
# Stabilization of the Final Sij_Boot matrix
Sij <- k10$Sij/k10$N
Sij <- Sij/max(Sij)

# Color palette
colC <- colorRampPalette(c("gold","darkgoldenrod2","seagreen","royalblue4","white","white","white","royalblue4","seagreen","darkgoldenrod2","gold"))
# Plot Stability matrix
par(mfrow=c(1,1))
try(corrplot(Sij,order="hclust",tl.col="black",method="color"
         ,addgrid.col=NA,col=colC(100),is.corr = FALSE,cl.lim = c(0,1)))

# ----------------------------------------------------------------------------------- #
#       Hierarchical Agglomerative Clustering
# ----------------------------------------------------------------------------------- #
# Final step: Hierarchical aglomerative clustering of the joint probability matrix
# http://scikit-learn.org/stable/modules/clustering.html

# Calculo de la matriz de distancias euclidenanas
vec <- rownames(Sij)
d <- dist(Sij,method = "euclidean")

# Método de Cluster
hc <- hclust(d,method = "ward.D",members = vec)
hc.dend <- as.dendrogram(hc)

# # Color the first sbranch in golden... the second sub-branch in green and the second sub-branch  in blue
# hc.dend[[1]] = dendrapply(hc.dend[[1]], colbranches, "darkgoldenrod2")
# hc.dend[[2]][[2]] = dendrapply(hc.dend[[2]][[2]], colbranches, "seagreen")
# hc.dend[[2]][[1]] = dendrapply(hc.dend[[2]][[1]], colbranches, "darkblue")
# 
# # Dendrogram
# par(mfrow=c(2,2))
# plot(hc.dend,xlab="",main="Dendrogram",col.main="black",cex.main=3,col="black",lwd=4,axes = FALSE,ylab="High")
# axis(2,col="black",lwd=2,at = seq(0,10,5),lab = seq(0,10,5),cex.axis=1.2,las=2,col.axis="black")
# 
# # Plots of each cluster
# # Selecciona solo 3 clusters
# clusterBest <- cutree(hc, 3)
# clust.ID<-as.data.frame(cbind(clusterBest,names(clusterBest)))
# colnames(clust.ID)<-c("clust","urm")
# Npsy.clust <- merge(Npsy.z,clust.ID,by="urm")
# x <- 1:10
# len <- length(Npsy.z)
# pts <- list(Npsy.z[,3:length(Npsy.z)])
# lab <- colnames(Npsy.z[,3:length(Npsy.z)])
# blank.plot(x,lab,"Cluster 1")
# plotLines(Npsy.clust[Npsy.clust$clust == "1",][3:12],"darkgoldenrod2",x)
# blank.plot(x,lab,"Cluster 2")
# plotLines(Npsy.clust[Npsy.clust$clust == "2",][3:12],"seagreen",x)
# blank.plot(x,lab,"Cluster 3")
# plotLines(Npsy.clust[Npsy.clust$clust == "3",][3:12],"darkblue",x)

# # ----------------------------------------------------------------------------------- #
# #               HAC Matrix of the clusterization
# # ----------------------------------------------------------------------------------- #
# par(mfrow=c(1,1))
# # Creates the Joint probability matrix & a vector that counts N-row selection by random resampling
# cM<-(matrix(data=0,nrow = n,ncol = n))
# colnames(cM)<-rownames(cM)<-names(clusterBest)
# 
# # Plot the HAC of the final cluster selection 
# cM<-cluster2mtx(clusterBest)
# colC <- colorRampPalette(c("white","white","white","white","darkgoldenrod2","seagreen","darkblue"))
# corrplot(cM,tl.col="black",method="color",addgrid.col="gray95",bg="black",is.corr = FALSE,
#          order="hclust",title = "",col=colC(200),cl.pos = 'n')
# 
