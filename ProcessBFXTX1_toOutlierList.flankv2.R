# KE Lotterhos
# Aug 10 2014
# Make Master List of bayenv2 results from Spruce



source("http://www.bioconductor.org/biocLite.R")
class(biocLite)
biocLite("limma")
library(limma)
install.packages("VennDiagram")
library(VennDiagram)

### Read in Sam's table ####
  ### Use these lines for spruce ####
  if(species=="SPRUCE"){
    wd <- "/data/seqcap/spruce/bwa_pseudo/final_tables/bayenv2Results10nonimputed_rm342479"
    setwd(wd)
    filepath<-"/data/seqcap/spruce/bwa_pseudo/"
    tablename <- "var_out_588ind_all_v2.summary.contig.ALL.annots.sorted.flank.pvalue_assoc"
    d1 <- read.table(paste(filepath,tablename,sep=""), header=TRUE, comment.char = "&")
      XTXcut <- read.table("Spruce_SeqCap.XTX.BFMAF.v2.XTXcutoffs", header=TRUE)
    head(XTXcut)
    BFcut <- read.table("Spruce_SeqCap.XTX.BFMAF.v2.BFcutoffs", header=TRUE)
    head(BFcut)
  }
  ### 

  ### Use these lines for pine ####
if(species=="PINE"){
    wd <- "/data/seqcap/pine/bwa_pseudo/final_tables/bayenv2Results10nonimputed_rm580plus"
    setwd(wd)
    filepath<-"/data/seqcap/pine/bwa_pseudo/"
    tablename <- "var_out_pine_all.summary.contig.ALL.annots.sorted.flank.pvalue_assoc"
    d1 <- read.table(paste(filepath,tablename,sep=""), header=TRUE, comment.char = "&")
      XTXcut <- read.table("Pine_SeqCap.XTX.BFMAF.v2.XTXcutoffs", header=TRUE)
    head(XTXcut)
    BFcut <- read.table("Pine_SeqCap.XTX.BFMAF.v2.BFcutoffs", header=TRUE)
    head(BFcut)
}
  ### 

### Load in environment cluster data and cutoffs for XTX/BF
  envi <- read.table("sprucePineEnvi.clusters", header=TRUE)
  head(envi)


  colnames(d1)[1:20]
  # Make a new variable with just nonsyn instead of specific AA change
    annot2<- as.character(d1$X.annotation)
    annot2[grep("nonsyn",d1$X.annotation)] <- "nonsyn"
    annot2[grep("non_tcontig", d1$X.annotation)] <- "intergenic"
    annot2[is.na(annot2)] <- "intergenic"
  table(annot2)
  sum(is.na(annot2))

  # Remove the SNPs that Sam said should be ignored
    remove <- c(which(annot2=="mismatch_altref"),
                which(annot2=="HAS_INDEL"),
                which(is.na(d1$XTX)))

  BFcols <- grep(".BF",colnames(d1))
  colnames(d1)[BFcols]
  r2 <- which(is.infinite(rowSums(d1[,BFcols])))
  length(r2)
  remove <- c(remove, r2)

    remove<-unique(remove)
    length(remove)
  # What is difference between FLANK and UTR?
    # Flanking sequences are upstream or downstream of the UTRs. UTRs are transcribed, while flanking regions are not. The flanking region presumably includes promoters.
  # "HAS_INDEL" does this mean the snp is an indel?  how is this coded?
    # these are errors somewhere in the annotation pipeline and should be ignored. Something went wrong (there are very few of them).
  # "mismatch_altref"
    # These have mismatching alleles between the transcriptome and the seqcap panel. Presume that they are errors.
  # "multi-allelic"
    # Triallelic or more (it wouldn't be possible to evaluate the alternative allele for the change amino acid)
  # "non_tcontig"
    # These are in intergenic regions, presumably (but could be intronic, unknown)
  # "unk_adj"
    # These are known to be adjacent to a gene, but unknown whether intron or flanking

  ### Write a new table with the filtered data ####
    d1<-data.frame(annot2, d1)
    d <- d1[-remove,]
   # write.table(d, "var_out_588ind_all_v2.summary.contig.ALL.annots.sorted.flank.pvalue_assoc.filtered", col.names=TRUE, row.names=FALSE)

  ### Calculate the percentages of each category for the whole dataset ####
    null.annot2 <- table(d$annot2)
    null.annot2.perc <- null.annot2/sum(null.annot2)
  
    getoutperc <- function(is.out){
      null.annot3 <- c(table(d$annot2[is.out]))
      as.numeric(null.annot3/sum(null.annot3))
    }

  ### To write to summary file
    #   dim(d.1) # loci in Sam's table
    #   dim(d2) # number loci after removing probable errors and NA for XTX
    #   null.annot2 # distribution of SNPs into categories



### Cutoffs for XTX 99% tail and BF for 99.9% tail ####

  #rm(d2, d1.1)


  #XTXcutoff1 <- XTXcut$XTX.99  # use this cutoff for transcripts > 1 hit
  #BFcutoff1 <- BFcut$BF.999    # use this cutoff for transcripts > 1 hit
  
  XTXcutoff1 <- sort(d$XTX)[round(length(d$XTX))*.99]
  BFcols <- grep(".BF",colnames(d))
  BF <- unlist(d[,BFcols])
  sum(is.infinite(BF))
  BFcutoff1 <- sort(BF)[round(length(BF)*.999)]

  XTX.logic.99 <- d$XTX > XTXcutoff1
  sum(XTX.logic.99)
  BF.logic.999 <- matrix(d[,BFcols]>BFcutoff1, nrow=nrow(d))
  sum(BF.logic.999)

  # Number of categories a SNP is an outlier in
  NumCat <- rowSums(cbind(XTX.logic.99, BF.logic.999))

  table(NumCat)
  d[which(NumCat>=21),BFcols]
  d[which(NumCat>=21),1:10]
  Is.Outlier.XTX99.BF999 <- NumCat>0
  Is.Outlier1 <- Is.Outlier.XTX99.BF999

### Outlier1: Explore enrichment for tcontigs with increasing # outliers ####
  get.trans.hits <- function(numhits){
      tc <- as.character(d$tcontig[Is.Outlier1])
      tc.table <- table(tc)
      trans.mult.hits <- names(tc.table)[which(tc.table>numhits)]
      length(tc.table); length(trans.mult.hits)
      which.mult <- (d$tcontig %in% trans.mult.hits)
      is.out <- (Is.Outlier1 & which.mult)
      list(is.out=is.out , numSNPs = sum(is.out), numtrans = length(trans.mult.hits))
  }
  #   Is.Outlier.trans.hits.gt5 <- get.trans.hits(5); lapply(Is.Outlier.trans.hits.gt5,sum)
    enrich2 <- NA
    colnames1 <- NA
    if (species=="SPRUCE"){thits <- c(0,1,5,10,15,20)}
    if (species=="PINE"){thits <- c(0,1,5,8, 10,13,14, 16,17)}
    for (transhits in thits){
        colnames1<-c(colnames1, transhits)
        print(transhits)
        Is.Outlier <- get.trans.hits(transhits)
        is.out <- Is.Outlier[[1]]
        stat<-Is.Outlier[2:3]
        change=c(round(-null.annot2.perc + getoutperc(is.out),5), stat)
        enrich2<-data.frame(enrich2,V=t(data.frame(change)))
      }
    colnames(enrich2)=colnames1
    enrich2
    write.table(enrich2, "EnrichmentByTransHits")
    trans.contig.ID <- as.numeric(as.factor(d$tcontig))
  ### LOOK AT THIS TCONTIG IN SPRUCE! ####
  if (species=="SPRUCE"){
  Is.Outlier.trans.hits.gt20 <-  get.trans.hits(20)
  d$tcontig[Is.Outlier.trans.hits.gt20[[1]]]
  table(d$annot2[Is.Outlier.trans.hits.gt20[[1]]])
    #this gene comp11186_c0_seq1 has 24 outlier SNPs, 13 of them are nonsyn and 7 in 5'UTR
  }
  ### LOOK AT THESE TCONTIGs IN PINE, each with 17 SNPs! ####
  if (species=="PINE"){
    Is.Outlier.trans.hits <- get.trans.hits(16)
    d$tcontig[Is.Outlier.trans.hits[[1]]]
    #comp12486_c0_seq1 15 in unk_flank, 2 intron
    #comp9540_c0_seq1  13 in 3' flank
    #comp16249_c0_seq1 11 unk_flank, 4 unk_ORF
   table(d$annot2[d$tcontig=="comp12486_c0_seq1" & Is.Outlier.trans.hits[[1]]])
   table(d$annot2[d$tcontig=="comp9540_c0_seq1" & Is.Outlier.trans.hits[[1]]])
   table(d$annot2[d$tcontig=="comp16249_c0_seq1" & Is.Outlier.trans.hits[[1]]])
  }

### Outlier2: Explore enrichment for outliers in increasing # categories ####
  Is.Outlier.numhits <- function(num.hits){
    outlier <- NumCat >= num.hits
    tc <- as.character(d$tcontig[outlier])
    tc.table <- table(tc)
    list(is.out = outlier, numSNPs = sum(outlier), numtrans = length(tc.table))
  }

  enrich3 <- NA
  for (num.hits in 1:max(NumCat)){
    outlist <- Is.Outlier.numhits(num.hits)
    is.out <- outlist[[1]]
      stat<-outlist[2:3]
      change=c(round(-null.annot2.perc + getoutperc(is.out),5), stat)
      enrich3<-data.frame(enrich3,V=t(data.frame(change)))
  }
  colnames(enrich3) <- c(NA, 1:max(NumCat))
  enrich3
  write.table(enrich3, "EnrichmentByNumCategories")

  if (species=="PINE"){
    Is.Outlier.trans.hits <- Is.Outlier.numhits(18)
    d$tcontig[Is.Outlier.trans.hits[[1]]]
    #5 are intergenic, 1 is in a 5'UTR
   table(d$annot2[Is.Outlier.trans.hits[[1]]])
  }

### Outlier3: Explore enrichment for outliers in more extreme parts of tail ####
  Is.Outlier.tail <- function(XTXcutoff2, BFcutoff2){
    XTX.logic <- d$XTX > XTXcutoff2
    BF.logic <- matrix(d[,BFcols]>BFcutoff2, nrow=nrow(d))
    outlier <- rowSums(cbind(XTX.logic, BF.logic))>0
    tc <- as.character(d$tcontig[outlier])
    tc.table <- table(tc)
    list(is.out = outlier, numSNPs = sum(outlier), numtrans = length(tc.table),
         num.XTX = sum(XTX.logic), num.BF = sum(rowSums(BF.logic)>0))
  }

  #d$tcontig[is.out2[[1]]]
  #table(d$annot2[is.out2[[1]]])
  #getoutperc(is.out2[[1]])-null.annot2.perc
 
  if (species=="SPRUCE"){
    xtx.tabs <- c(XTXcutoff1, 1000,1600,2000,2200, 2500)
    bf.tabs <- c(BFcutoff1, 30, 45, 70, 90, 110)
  }
  if(species=="PINE"){
    xtx.tabs <- c(XTXcutoff1, 5000, 7000, 8700, 10000,12500,15000)
    bf.tabs <- c(BFcutoff1, 300, 400, 500, 600, 700)
  }
 enrich <- NA
  colnames1 <- NA
  for (XTX in xtx.tabs){
    for (BF in bf.tabs){
      colnames1<-c(colnames1, paste(XTX,BF,sep="_"))
      print(paste(XTX,BF,sep="_"))
      Is.Outlier <- Is.Outlier.tail(XTX, BF)
      is.out <- Is.Outlier[[1]]
      stat<-Is.Outlier[2:5]
      change=c(round(-null.annot2.perc + getoutperc(is.out),5), stat)
      enrich<-data.frame(enrich,V=t(data.frame(change)))
    }
  }
  colnames(enrich)=colnames1
  enrich
  write.table(enrich, "EnrichmentByCutoff")


### Suggested cutoffs for both species ####
  # In 99% tail, a transcript with > 5 hits
    Is.Outlier.trans.hits.gt5 <- get.trans.hits(5)
  # In 99% tail, a SNP outlier in > 10 factors
    Is.Outlier.BF.hits.gt10 <- Is.Outlier.numhits(10)


  if (species=="SPRUCE"){
    XTX.cut <- 1600
      BF.cut <- 45
    Is.Outlier.XTX.BF <- Is.Outlier.tail( XTX.cut,BF.cut)
  }
  if (species=="PINE"){
        XTX.cut <- 8700
      BF.cut <- 500
    Is.Outlier.XTX.BF <- Is.Outlier.tail( XTX.cut,BF.cut)
  }

  Is.Outlier.Overall <- Is.Outlier.trans.hits.gt5[[1]] | Is.Outlier.BF.hits.gt10[[1]] |
             Is.Outlier.XTX.BF[[1]]
  (round(-null.annot2.perc + getoutperc(Is.Outlier.Overall),5))
  


   GetSnpsTranscripts <- function(is.out){
      c(sum(is.out), length(table(as.character(d$tcontig[is.out]))))
   }
   GetSnpsTranscripts(Is.Outlier.Overall)
  

  vennCounts(cbind(Is.Outlier.trans.hits.gt5[[1]], Is.Outlier.BF.hits.gt10[[1]],
             Is.Outlier.XTX.BF[[1]]))

### Outliers by Environmental Group
    BF.logic <- matrix(d[,BFcols]>BF.cut, nrow=nrow(d))
    Is.Outlier.XTX = d$XTX > XTX.cut
  #XTX.logic.99 
  #BF.logic.999
#####
    colnames(BF.logic) <- paste(colnames(d[BFcols]), ".Is.Outlier", sep="")
    colnames(BF.logic.999) <- paste(colnames(d[BFcols]), ".Is.999Outlier")

    Is.Outlier.G1Temp <- (rowSums(BF.logic[,envi$G3==1])>0) 
    Is.Outlier.G1Temp.999 <- (rowSums(BF.logic.999[,envi$G3==1])>0)
   #  (round(-null.annot2.perc + getoutperc(Is.Outlier.G1Temp),5))
         GetSnpsTranscripts(Is.Outlier.G1Temp)
          GetSnpsTranscripts(Is.Outlier.G1Temp.999)

    Is.Outlier.G2Precip <- (rowSums(BF.logic[,envi$G3==2])>0)
    Is.Outlier.G2Precip.999 <- (rowSums(BF.logic.999[,envi$G3==2])>0)
   #  (round(-null.annot2.perc + getoutperc(Is.Outlier.G2Precip),5))
         GetSnpsTranscripts(Is.Outlier.G2Precip)
         GetSnpsTranscripts(Is.Outlier.G2Precip.999)

    Is.Outlier.G3Frost <- (rowSums(BF.logic[,envi$G3==3])>0)
    Is.Outlier.G3Frost.999 <- (rowSums(BF.logic.999[,envi$G3==3])>0)
    #  (round(-null.annot2.perc + getoutperc(Is.Outlier.G2Precip),5))
          GetSnpsTranscripts(Is.Outlier.G3Frost)
          GetSnpsTranscripts(Is.Outlier.G3Frost.999)

  #sum(Is.Outlier.G1Temp | Is.Outlier.G2Precip | Is.Outlier.G3Frost)
  #sum(rowSums(BF.logic)>0)
  #vennCounts(cbind(Is.Outlier.G1Temp, Is.Outlier.G2Precip,
  #           Is.Outlier.G3Frost)
             
### Write Tables ####

  d.new <- data.frame(d,
                      NumCat.XTX.BF = NumCat,
                      Is.Outlier.XTX99.BF999,
                      Is.Outlier.trans.hits.gt5=Is.Outlier.trans.hits.gt5[[1]],
                      Is.Outlier.BF.hits.gt10=Is.Outlier.BF.hits.gt10[[1]],
                      Is.Outlier.XTX.BF=Is.Outlier.XTX.BF[[1]],
                      Is.Outlier.Overall,
                      Is.Outlier.G1Temp,Is.Outlier.G1Temp.999,
                      Is.Outlier.G2Precip,Is.Outlier.G2Precip.999,
                      Is.Outlier.G3Frost,Is.Outlier.G3Frost.999,
                      Is.Outlier.XTX, XTX.logic.99,
                      BF.logic, BF.logic.999)
  write.table(d.new[Is.Outlier.XTX99.BF999,],paste(tablename, ".filtered.outliers", sep=""),  
              col.names=TRUE, row.names=FALSE)     
write.table(d.new, paste(tablename, ".filtered", sep=""), 
              col.names=TRUE, row.names=FALSE)
       

   Top10 <- Is.Outlier.trans.hits.gt5[[1]] & Is.Outlier.BF.hits.gt10[[1]] &
              Is.Outlier.XTX.BF[[1]]
   sum(Top10)
   c(table(d$annot2[Top10]), sum(is.na(d$annot2[Top10])))
#   round(-null.annot2.perc + getoutperc(Top10),5)
# 
   Top100 <- (Is.Outlier.BF.hits.gt10[[1]] & Is.Outlier.XTX.BF[[1]]) | 
             (Is.Outlier.trans.hits.gt5[[1]] & Is.Outlier.XTX.BF[[1]]) |
             (Is.Outlier.BF.hits.gt10[[1]] & Is.Outlier.trans.hits.gt5[[1]])
   sum(Top100)
   c(table(d$annot2[Top100]))
   round(-null.annot2.perc + getoutperc(Top100),5)
#   round(-null.annot2.noNA.perc + getoutperc.noNA(Top100),5)





#   ### For output to report
#     sum(XTX.logic2) # number SNPs
#     length(unique(d$tcontig[XTX.logic2]))  # number contigs
#     length(unique(d$tcontig))
#     BF.999.logic.overall <- (rowSums(BF.logic2)>0) 
#     sum(BF.999.logic.overall)  # number SNPs across all environments
#     length(unique(d$tcontig[BF.999.logic.overall])) # number contigs across all envis
#     colnames(BF.999.logic) <- paste(colnames(d)[BFcols], ".gt.999", sep="")
#     colSums(BF.999.logic) # number in each environment  

  # Get logic for group 1 or 2, for 2cluster
  # envi$G2 is group of BF into 2 clusters

#   BF.999.twoG1 <- (rowSums(BF.999.logic[,envi$G2==1])>0) 
#   BF.999.twoG2 <- (rowSums(BF.999.logic[,envi$G2==2])>0) 
#   
#   # envi$G6 is group of BF into 6 clusters
#   BF.999.sixG1 <- (rowSums(BF.999.logic[,envi$G6==1])>0) 
#   BF.999.sixG2 <- (rowSums(BF.999.logic[,envi$G6==2])>0) 
#   BF.999.sixG3 <- (rowSums(BF.999.logic[,envi$G6==3])>0) 
#   BF.999.sixG4 <- (rowSums(BF.999.logic[,envi$G6==4])>0) 
#   BF.999.sixG5 <- (rowSums(BF.999.logic[,envi$G6==5])>0) 
#   BF.999.sixG6 <- (rowSums(BF.999.logic[,envi$G6==6])>0) 


  
 




