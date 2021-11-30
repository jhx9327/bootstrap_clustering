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
home <- "C:\\1-hongxiu\\codes"

# Load data
setwd(paste0(home,"\\2020_cognition_connectomics_TLE\\code\\BASC\\"))
source("functions.R")
source("BASC_func.R")
Npsy <- read.csv("Used_asd_control.csv")

# Erases the rows with NA values
#Npsy <- Npsy[complete.cases(Npsy),vars] 

# Selects only those variables of interest
Npsy.z <- Npsy[,c("site","FILE_ID","DX_GROUP","AGE_AT_SCAN","ADOS_TOTAL","ADOS_COMM","ADOS_SOCIAL","ADOS_STEREO_BEHAV")]

# Create a Z-score based on controls for each cognitive variable
controls <- which(Npsy.z$DX_GROUP == 2)
asd <- which(Npsy.z$DX_GROUP == 1)
# Npsy.z[,4:8] <- apply(Npsy.z[,4:8], 2, function(x) Zscore(x[asd],x)) # 2 means by row, 1 means by column
Npsy.z[asd,5:8]<- apply(Npsy.z[asd,5:8], 2, function(x) Zscore_hx(x[asd])) # 2 means by row, 1 means by column

#### Clusterizacion Jerarquica con Matriz de distancias euclidianas ####
# Matrix de valores normalizados
# Npsy.matrix <- as.matrix(Npsy.z[Npsy.z$group != "ctrl",][3:12])


# ----------------------------------------------------------------------------------- #
#   Visualization of z-scored dataset
# ----------------------------------------------------------------------------------- #
# https://www.rdocumentation.org/packages/gplots/versions/3.0.1/topics/heatmap.2
Npsy.matrix <- as.matrix(Npsy.z[asd,5:8])
n <- dim(Npsy.matrix)[1]
rownames(Npsy.matrix)<-sprintf("ASD %02d", c(1:n))
colC <- colorRampPalette(c("darkblue","royalblue4","seagreen","darkgoldenrod2","gold"))
heatmap.2(Npsy.matrix,col=colC(256),
          scale = "none",key=TRUE,key.title = "",denscol = "black",tracecol = NA,dendrogram = "none",main = "age and syptom scores",Rowv=FALSE,Colv=FALSE)


# ----------------------------------------------------------------------------------- #
#   Visualization of the Cluster Selection Procedure
# ----------------------------------------------------------------------------------- #
# Row ID with URM
# rownames(Npsy.matrix)<-Npsy.z[Npsy.z$group != "ctrl",][[1]]

####### Total within-clusters sum of squares (wss) for k-means 
asd_msn_corrected <- as.matrix(read.csv('meanMS_regional_corrected_asd.csv', header = FALSE))
Npsy.matrix <- asd_msn_corrected
# par(mfrow=c(1,3))
k.max <- 20#dim(Npsy.matrix)[2] 
wss <- sapply(1:k.max,
              function(k){kmeans(Npsy.matrix, k, nstart=50,iter.max = 15 )$tot.withinss})
plot(1:k.max, wss,las=2,lwd=3,axes=FALSE,cex.lab=2,cex.main=2.25,
     type="b", pch = 19, frame = FALSE,
     xlab="Number of Clusters",
     ylab="",main="Within-Clusters Sum of Squares") # WCSS
axis(1,at=seq(2,k.max,2),labels = TRUE,cex.axis=1.5,lwd=1.5)
axis(2,at=seq(100,500,200),labels = TRUE,las=2,cex.axis=1.5,lwd=1.5)

#######  Clusterization with ward.D2 and euclidean distance with Best number of clusters
clus <- MyNbClust(Npsy.matrix, distance = "euclidean", min.nc=2, max.nc=8, 
             method = "ward.D2", index = "all",plotetc = FALSE) # min.nc=2, 
clusterBest <- clus$Best.partition
k <- max(clusterBest)
abline(v=k,lty=2,lwd=2,col="red")

# Best Cluster partition selection
barplot(table(clus$Best.nc[1,]),border = "navy",col=addTrans("navy",65),space=0,cex.lab=2,
        xlab="Number of clusters", ylab="",axes=F,xaxt='n',main="Number of clusters selected by 26 criteria")
axis(1,at=seq(0.5,7.5,1),labels = c(0:4,6:8),cex.axis=1.5,lwd=1.5,pos = -0.25)
axis(2,at=seq(0,8,2),labels = TRUE,cex.axis=1.5,lwd=1.5)
mtext("Frequency of selection",2,line = 2.5,at = 4,cex = 1)
points(3.5,8,pch=18,col="black",cex=4)
points(3.5,8,pch=18,col="red",cex=3)

# Multi Dimensional Scaling Plot
col1 <- colorRampPalette(c("firebrick1","forestgreen","dodgerblue3"))
d <- dist(Npsy.matrix,method = "euclidean")
distances <- as.matrix(d)
mds = isoMDS(d)
plot(mds$points,las=1,bty='n',pty='l',pch = 21, cex = 5, col=col1(k)[clusterBest],bg = adjustcolor(col1(k)[clusterBest], alpha = 0.2), 
     xlab = "X", ylab = "Y",xlim=c(-12,10),ylim=c(-5,4)) 
text(mds$points, labels = rownames(mds$points), cex = 0.7,col=col1(k)[clusterBest])


# ----------------------------------------------------------------------------------- #
#       A BINARY SIMILARITY MATRIX
# ----------------------------------------------------------------------------------- #
clusterBest<-clus$Best.partition
s<-names(clusterBest)
names(clusterBest)<-s<-sapply(strsplit(s, split='.', fixed=TRUE), function(x) (x[1]))
# Stability matrix from original sample
S0 <- (matrix(data=0,nrow = n,ncol = n)) # makes and empty matrix
colnames(S0) <- rownames(S0) <- rownames(Npsy.matrix)
S.orig<-mtx.jointProb(clusterBest,S0)
# Plot of the similarity matrix of the orginal data
par(mfrow=c(1,1))
corrplot(S.orig,order="alphabet",tl.col="black",method="color",addgrid.col="gray95",bg="black",is.corr = FALSE)
# Removes unnecesary variables
rm(clus,clusterBest,d,k,k.max,mds,wss,s,S0,S.orig)


# ----------------------------------------------------------------------------------- #
#       BASC - Boostrap Analysis of Stable Clusters
# ----------------------------------------------------------------------------------- #
k10 <- BASC(Matrix = Npsy.matrix, rows.ID = rownames(Npsy.matrix), boots = 10000)

# Saves the result of the 10,000 bootstrap to a .RData file
#save(k10,file = paste0(home,"/git_here/micasoft/sandbox/raul/BASC/10k_34subjs.RData"))

#Loads the data
load(paste0(home,"/git_here/micasoft/sandbox/raul/BASC/10k_34subjs_new.RData"))

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

# Color the first sbranch in golden... the second sub-branch in green and the second sub-branch  in blue
hc.dend[[1]] = dendrapply(hc.dend[[1]], colbranches, "darkgoldenrod2")
hc.dend[[2]][[2]] = dendrapply(hc.dend[[2]][[2]], colbranches, "seagreen")
hc.dend[[2]][[1]] = dendrapply(hc.dend[[2]][[1]], colbranches, "darkblue")

# Dendrogram
par(mfrow=c(2,2))
plot(hc.dend,xlab="",main="Dendrogram",col.main="black",cex.main=3,col="black",lwd=4,axes = FALSE,ylab="High")
axis(2,col="black",lwd=2,at = seq(0,10,5),lab = seq(0,10,5),cex.axis=1.2,las=2,col.axis="black")

# Plots of each cluster
# Selecciona solo 3 clusters
clusterBest <- cutree(hc, 3)
clust.ID<-as.data.frame(cbind(clusterBest,names(clusterBest)))
colnames(clust.ID)<-c("clust","urm")
Npsy.clust <- merge(Npsy.z,clust.ID,by="urm")
x <- 1:10
len <- length(Npsy.z)
pts <- list(Npsy.z[,3:length(Npsy.z)])
lab <- colnames(Npsy.z[,3:length(Npsy.z)])
blank.plot(x,lab,"Cluster 1")
plotLines(Npsy.clust[Npsy.clust$clust == "1",][3:12],"darkgoldenrod2",x)
blank.plot(x,lab,"Cluster 2")
plotLines(Npsy.clust[Npsy.clust$clust == "2",][3:12],"seagreen",x)
blank.plot(x,lab,"Cluster 3")
plotLines(Npsy.clust[Npsy.clust$clust == "3",][3:12],"darkblue",x)

# ----------------------------------------------------------------------------------- #
#               HAC Matrix of the clusterization
# ----------------------------------------------------------------------------------- #
par(mfrow=c(1,1))
# Creates the Joint probability matrix & a vector that counts N-row selection by random resampling
cM<-(matrix(data=0,nrow = n,ncol = n))
colnames(cM)<-rownames(cM)<-names(clusterBest)

# Plot the HAC of the final cluster selection 
cM<-cluster2mtx(clusterBest)
colC <- colorRampPalette(c("white","white","white","white","darkgoldenrod2","seagreen","darkblue"))
corrplot(cM,tl.col="black",method="color",addgrid.col="gray95",bg="black",is.corr = FALSE,
         order="hclust",title = "",col=colC(200),cl.pos = 'n')

