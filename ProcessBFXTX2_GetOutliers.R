species<-"SPRUCE"

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

### Get top transcript hits
  get.top.trans <- function(num){
    cont <- as.character(d$tcontig[d$Is.Outlier.trans.hits.gt5])
    sort.cont <- sort(table(cont),decreasing = TRUE)
    ind<- which(sort.cont>=sort.cont[num]) #get ties
    outlier.names <- names(sort.cont)[ind]
    d[which(d$tcontig %in% outlier.names),]
  }
  
  top.trans.hits <- get.top.trans(1)
  top.trans.hits[,1:10]

### Get SNPs with multiple environment hits
  get.top.mult.envi <- function(num){
    # works slightly differently than last function
    cont <- d$NumCat.XTX.BF[d$Is.Outlier.BF.hits.gt10]
    sort.cont <- rev(table(cont))
    ind<- 1:num #get ties
    d[which(d$NumCat.XTX.BF %in% names(sort.cont)[ind]),]
  }
  
  top.mult.envi.hits <- get.top.mult.envi(2)
  top.mult.envi.hits[,1:10]

### Get top XTX
  get.top.XTX <- function(num){
    d[order(d$XTX,decreasing = TRUE)[1:num],]
  }
  top.XTX.hits <- get.top.XTX(5)
  top.XTX.hits[,1:10]
  top.XTX.hits$XTX

### Get top BF in a specific environment
 get.top.BF.envi <- function(num, envi.name){
    column <- grep(envi.name, colnames(d))[1]
    d[order(d[,column],decreasing = TRUE)[1:num],]
  }
  elev.hits <- get.top.BF.envi(3, "ELEVATION")
  elev.hits[,1:10]
  if(species=="SPRUCE"){elev.hits$ELEVATION.BF}
  if(species=="PINE"){elev.hits$ELEVATION.BF_rank}

### Get top BF across all environments
    BFcols <- grep(".BF", colnames(d))[1:22]
    max.BF <- apply(d[,BFcols],1,max)
    d <- data.frame(d, max.BF) 
    get.top.BF.overall <- function(num){
      d[order(max.BF, decreasing=TRUE)[1:num],]
    }
    # Keep in mind that some environments have very long BF tails, so this 
    # will bias SNPs that are outliers in those environments
    BF.hits <- get.top.BF.overall(5)
    BF.hits[,1:10]
    BF.hits$max.BF
