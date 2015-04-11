# CalcEnviClusters.R
# KE Lotterhos
# Oct 22, 2014

# Takes a matrix of environmental variables across populations as input
# Does a cluster analysis
# Assigns environmental variables to a cluster

setwd("/Users/katie/Desktop/AdaptreeData/src_pine/src-remove580-582,584-586")
envimat <- "var_out_pine_all_COMBINED.table.contig_flt10.bayenv.envi2"
envinames <- "enviNamesAllAnalyses.txt"

  envi <- read.table(envimat)
  rowenvi <- read.table(envinames)
  dim(envi)
  rownames(envi) <- rowenvi$V1
  e <- scale(envi)
  # K-Means Cluster Analysis
  #fit3 <- kmeans(e, 3) # 3 cluster solution
  # get cluster means 
  #clusterMeans <- t(aggregate(e,by=list(fit3$cluster),FUN=mean))

  # Ward Hierarchical Clustering
  makeplot<- function(d, labs, k1){
    fit <- hclust(d, method="ward.D") 
    pdf(paste(envimat,".",k1, "clusters.pdf", sep=""), width=6, height=8)
      plot(fit, labels = labs, main="", xlab="") # display dendogram
      groups <- cutree(fit, k=k1) # cut tree into x clusters
      # draw dendogram with red borders around the x clusters 
      rect.hclust(fit, k=k1, border="red")
    dev.off()
    return(as.numeric(groups))
  }
  par(mfrow=c(1,1), mar=c(1,4,1,1))
  d <- dist(e, method = "euclidean")
  G6 <- makeplot(d, rowenvi$V2, 6)
  G2 <- makeplot(d, rowenvi$V2, 2)
  # append cluster assignment
  mydata <- data.frame(enviAbb = rowenvi$V1, enviDesc = rowenvi$V2, G2, G6)
  mydata$G6.desc <- NA
  mydata$G6.desc[G6==1] <- "Lat and freezing"
  mydata$G6.desc[G6==2] <- "Long, summer precip"
  mydata$G6.desc[G6==3] <- "Elev, prec, frost-free"
  mydata$G6.desc[G6==4] <- "frost-free"
  mydata$G6.desc[G6==5] <- "Evap, max temp"
  mydata$G6.desc[G6==6] <- "Heat:moisture"
  
write.table(mydata, paste(envimat, ".clusters", sep=""), row.names=FALSE, col.names=TRUE)
