
# source("http://www.bioconductor.org/biocLite.R")
# class(biocLite)
# biocLite("limma")
library(limma)
#install.packages("VennDiagram")
library(VennDiagram)

if (species=="SPRUCE"){
  setwd("/Users/katie/Desktop/AdaptreeData/results_spruce")
  d<- read.table("var_out_588ind_all_v2.summary.contig.ALL.annots.sorted.flank.pvalue_assoc.filtered.outliers", header=TRUE)
}
if (species=="PINE"){
   setwd("/Users/katie/Desktop/AdaptreeData/results_pine")
   d <- read.table("var_out_pine_all.summary.contig.ALL.annots.sorted.flank.pvalue_assoc.filtered.outliers", header=TRUE)
}


e<-read.table("~/Desktop/AdaptreeData/src_PineSpruceFunctions/envi clusters/sprucePineEnvi.clusters", header=TRUE)
colnames(d)
numtran <- length(levels(d$tcontig))

### Histogram of SNPs per transcript
  pdf("Hist_SnpsPerTranscContig.pdf", width=10, height=4)
    yo<-table(as.character(d$tcontig[d$Is.Outlier.XTX99.BF999]))
    par(mar=c(4,4,2,1))
    plot(yo[order(yo,decreasing = TRUE)], xaxt = "n", type="h",
         xlab = "Transcriptomic contig", ylab="Number of outlier SNPs in 99.9% tail",
         main=species)
    axis(1, at=c(0,numtran/2,numtran))
  dev.off()

### Venn counts for 3-cluster environment ####  
  test1 <- cbind(d$Is.Outlier.XTX, 
                 d$Is.Outlier.G1Temp, 
                 d$Is.Outlier.G2Precip,
                 d$Is.Outlier.G3Frost)
  pdf("VennDiagramXTX.Envigroups.SNPs.pdf", width=12, height=12)
    vennDiagram(vennCounts(test1), names = c("XTX", "Temp",
                             "Precipitation",
                             "Frost"), 
                main=paste(species, "SNPs"), cex=2)
# circle.col=c("orange", "orange","blue","red","black"), cex=3)
  dev.off()

  pdf("VennDiagramEnvigroups.SNPs.pdf", width=12, height=12)
    vennDiagram(vennCounts(cbind( 
                 d$Is.Outlier.G1Temp, 
                 d$Is.Outlier.G2Precip,
                 d$Is.Outlier.G3Frost)
      ), names = c("Temp",
                             "Precipitation",
                             "Frost"), 
                main=paste(species, "SNPs"), cex=2)
# circle.col=c("orange", "orange","blue","red","black"), cex=3)
  dev.off()

  pdf("VennDiagramEnvigroups.SNPs999.pdf", width=12, height=12)
    vennDiagram(vennCounts(cbind( 
                 d$Is.Outlier.G1Temp.999, 
                 d$Is.Outlier.G2Precip.999,
                 d$Is.Outlier.G3Frost.999)
      ), names = c("Temp",
                             "Precipitation",
                             "Frost"), 
                main=paste(species, "SNPs"), cex=2)
# circle.col=c("orange", "orange","blue","red","black"), cex=3)
  dev.off()

### Venn for transcript contigs: BF groups, XTX ####
  BF1 <- table(as.character(d$tcontig),d$Is.Outlier.G1Temp) 
  BF2 <- table(as.character(d$tcontig),d$Is.Outlier.G2Precip)
  BF3 <- table(as.character(d$tcontig), d$Is.Outlier.G3Frost)
  XTX <- table(as.character(d$tcontig),d$Is.Outlier.XTX)
  
  test2<- cbind(XTX[,2]>0,BF1[,2]>0,BF2[,2]>0, BF3[,2]>0)

  pdf("VennDiagramXTX.Envigroups.transcripts.pdf", width=12, height=12)
    vennDiagram(vennCounts(test2), names = c("XTX", "Temp",
                             "Precipitation",
                             "Frost"),  
                main=paste(species, "transcripts"), cex=2)
# circle.col=c("orange", "orange","blue","red","black"), cex=3)
  dev.off()

  pdf("VennDiagramEnvigroups.transcripts.pdf", width=12, height=12)
    vennDiagram(vennCounts(cbind(BF1[,2]>0,BF2[,2]>0, BF3[,2]>0))
      , names = c("Temp",
                             "Precipitation",
                             "Frost"),  
                main=paste(species, "transcripts"), cex=2)
# circle.col=c("orange", "orange","blue","red","black"), cex=3)
  dev.off()

  BF1.999 <- table(as.character(d$tcontig),d$Is.Outlier.G1Temp.999) 
  BF2.999 <- table(as.character(d$tcontig),d$Is.Outlier.G2Precip.999)
  BF3.999 <- table(as.character(d$tcontig), d$Is.Outlier.G3Frost.999)
  
  test3<- cbind(BF1.999[,2]>0,BF2.999[,2]>0, BF3.999[,2]>0)
  pdf("VennDiagramEnvigroups.transcripts999.pdf", width=12, height=12)
    vennDiagram(vennCounts(test3)
      , names = c("Temp",
                             "Precipitation",
                             "Frost"),  
                main=paste(species, "transcripts"), cex=2)
# circle.col=c("orange", "orange","blue","red","black"), cex=3)
  dev.off()

### what is overlap of BF for correlated envi variables? ####
    colnames(d)
    if (species=="SPRUCE"){BFcols <- grep(".BF.Is.Outlier", colnames(d))}
    if (species=="PINE"){BFcols <- grep(".BF_rank.Is.Outlier", colnames(d))}
    GroupID <- c("Temp", "Precipitation", "Frost")
    for (i in 1:3){
      ind<-e$G3==i
      z<- cbind(d[,BFcols[ind]])
      colnames(z)<- c(as.character(e$enviDesc[ind]))
      vcz <- vennCounts(z)
      numcat <- rowSums(vcz[,1:sum(ind)])
      vcz.ord <- order(numcat)
      write.table(vcz[vcz.ord,], 
                  paste("VennCountsGroup", GroupID[i],sep="")
      )
      print(c(sum(vcz[numcat==1,sum(ind)+1]),
              sum(vcz[numcat>1,sum(ind)+1]),
              sum(vcz[numcat==sum(ind),sum(ind)+1])
      ))
    }
#   pdf(paste("VennDiagramXTX.Group",i,".pdf",sep=""), width=12, height=12)
#   par(mar=c(1,1,6,1))
#   vennDiagram(vennCounts(z), cex=2, main=paste("Group", i))
#   dev.off()
#}
  if (species=="SPRUCE"){BFcols999 <- grep(".BF..Is.999Outlier", colnames(d))}
    if (species=="PINE"){BFcols999 <- grep(".BF_rank..Is.999Outlier", colnames(d))}
    GroupID <- c("Temp", "Precipitation", "Frost")
    for (i in 1:3){
      ind<-e$G3==i
      z<- cbind(d[,BFcols999[ind]])
      colnames(z)<- c(as.character(e$enviDesc[ind]))
      vcz <- vennCounts(z)
      numcat <- rowSums(vcz[,1:sum(ind)])
      vcz.ord <- order(numcat)
      write.table(vcz[vcz.ord,], 
                  paste("VennCounts999Group", GroupID[i],sep="")
      )
      print(c(sum(vcz[numcat==1,sum(ind)+1]),
              sum(vcz[numcat>1,sum(ind)+1]),
              sum(vcz[numcat==sum(ind),sum(ind)+1])
      ))
    }


### Lattice Plot in rho ####
  library(lattice)
  envi <- read.table("seqcap.bayenv.BF.envi2")
  dim(envi)
  envi <- t(envi)
  colnames(envi)<-e$enviAbb

  #rhocols <- grep(".rho",colnames(d))
    if (species=="SPRUCE"){BFcols2 <- grep(".BF", colnames(d))[1:22]}
    if (species=="PINE"){BFcols2 <- grep(".BF_rank", colnames(d))[1:22]}
  colnames(d)[BFcols2]
  e3names<- e$G3.desc
  ord<- order(e$G3)

 png(paste("ScatterAll.envi.png",sep=""), width=14, height=14, units = "in", res=300)
    splom(envi[,ord])
  dev.off()

#  png(paste("ScatterAll.rho.png",sep=""), width=14, height=14, units = "in", res=300)
#     splom(d[,rhocols[ord]])
#   dev.off()

 png(paste("ScatterAll.BF.png",sep=""), width=14, height=14, units = "in", res=300)
    splom(d[d$Is.Outlier.XTX.BF, ord+(BFcols2[1]-1)])
  dev.off()

 pdf(paste("CorrEnvi.BF.pdf",sep=""), width=4, height=4)
  # remove infinites 
  plot(cor(envi), 
       cor(d[d$Is.Outlier.XTX.BF,BFcols2]),
       xlab= "Correlation between 2 environmental variables",
       ylab= "Correlation between Bayes Factors of outliers")
  dev.off()

 pdf(paste("CorrEnvi.BF999.pdf",sep=""), width=4, height=4)
  # remove infinites 
  plot(cor(envi), 
       cor(d[d$Is.Outlier.XTX99.BF999,BFcols999]),
       xlab= "Correlation between 2 environmental variables",
       ylab= "Correlation between Bayes Factors of 99.9% outliers")
  dev.off()



### Ex of SNP adapting to multiple envis ####
  if (species=="SPRUCE"){
      envi2 <- data.frame(envi)
      pdf("ExampleOfSNPadaptingToMultipleEnvi.pdf", width=8, height=4)
      par(mfrow=c(1,2), mar=c(4,4,1,1))
      plot(envi2$MAP, envi2$FFP, bty="n",
           xlab="Mean Annual Precipitation",
           ylab="Frost-Free Period")
      abline(lm(envi2$FFP~envi2$MAP))
      plot(d$MAP.BF[d$Is.Outlier.XTX.BF],
           d$FFP.BF[d$Is.Outlier.XTX.BF],
           xlab= "Bayes Factor for MAP",
           ylab= "Bayes Factor for FFP",
           bty="n")
      dev.off()
    
    ### Ex of SNP tradeoff in correlated envis ####
      envi2 <- data.frame(envi)
      pdf("ExampleOfSNPtrade-off-CorrelatedEnvi.pdf", width=8, height=4)
      par(mfrow=c(1,2), mar=c(4,4,1,1))
      plot(envi2$MAP, envi2$MSP, bty="n",
           xlab="Mean Annual Precipitation",
           ylab="Mean Summer Precipitation")
      cor(envi2$MAP, envi2$MSP)
      abline(lm(envi2$MSP~envi2$MAP))
      cor(d$MAP.BF[d$Is.Outlier.XTX.BF],
           d$MSP.BF[d$Is.Outlier.XTX.BF])
      plot(d$MAP.BF[d$Is.Outlier.XTX.BF],
           d$MSP.BF[d$Is.Outlier.XTX.BF],
           xlab= "Bayes Factor for MAP",
           ylab= "Bayes Factor for MSP",
           bty="n")
      dev.off()
}

#### Venn diagrams for other cutoff criteria ####
colnames(d)
  pdf("VennDiagramOutliersSNPs.pdf", width=12, height=12)
    a<- vennCounts(cbind(d$Is.Outlier.trans.hits.gt5,
                         d$Is.Outlier.BF.hits.gt10,
                         d$Is.Outlier.XTX.BF)
                   )
    vennDiagram(a, names = c("In transcript with\n> 5 outliers", 
                             "Outlier in >10\nenvironments",
                             "Extreme tail"),
      main=species, cex=2)
# circle.col=c("orange", "orange","blue","red","black"), cex=3)
  dev.off()

  sum(d$Is.Outlier.trans.hits.gt5)
  sum(d$Is.Outlier.BF.hits.gt10)
  sum(d$Is.Outlier.XTX.BF)
  sum(d$Is.Outlier.XTX99.BF999)
# i=6 # for some reason a for loop or function doesn't work here
#   ind<-e$G6==i
#   print(i)
#   z.rho<- d[,rhocols[ind]]
#   z.bf <- d[,BFcols2[ind]]
#   png(paste("ScatterG",i,".rho.png",sep=""), width=10, height=10, units = "in", res=300)
#     splom(z.rho)
#   dev.off()
#   png(paste("ScatterG",i,"BF.png",sep=""), width=10, height=10, units = "in", res=300)
#     splom(z.bf)
#   dev.off()
#   png(paste("ScatterG",i,"envi.png",sep=""), width=5, height=5, units = "in", res=300)
#     splom(envi[,ind])
#   dev.off()

