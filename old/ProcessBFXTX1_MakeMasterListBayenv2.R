# KE Lotterhos
# Aug 10 2014
# Make Master List of bayenv2 results from Spruce or Pine


#ArrangeMasterList <- function(wkdir, locinames, XTXBF, MAF, envi){

	setwd(wkdir)

# Read in list of locus names
  loci.names <- read.table(locinames, 
                           colClasses = c("numeric","character"),
                          col.names=c("Num", "Name"))
  outfile <- paste(XTX.BF, ".v2")
	#str(loci.names)

# Read in list of environments  
  envi.df <- read.table(envi, header=FALSE)
	envi.names <- envi.df$V1
  xtx.bf <- read.table(XTX.BF,
                      colClasses=c("character", "numeric", "character", rep("numeric", length(envi.names)*3)))

# Check to make sure all loci names are equal
	writeLines("Is the list of loci names in XTX and BF identical?")
	print(identical(xtx.bf[,1], xtx.bf[,3]))
	
  envi.names2 <- (matrix(c(paste(envi.names,".BF",sep=""), 
                          paste(envi.names,".rho1",sep=""),  
                          paste(envi.names,".rho2",sep="")
                          ), ncol=3))
  envi.names3 <- as.vector(t(envi.names2))
  colnames(xtx.bf) <- c("Name1", "XTX", "Name2", envi.names3)
  #  xtx.bf$XTX<-as.numeric(as.character(xtx.bf$XTX))
	#str(xtx.bf$XTX)
	#str(xtx.bf[,1:5])

  repeated.loci <- loci.names$Name[duplicated(loci.names$Name)]
  writeLines("If the loci names has a repeated locus, print here:")
	print(repeated.loci)
	
  
  # remove non-unique lines
  if (length(repeated.loci)>0 ){
    which.in.xtx <-  xtx.bf$Name1 %in% repeated.loci
    if (length(which.in.xtx)>0){
       writeLines("Removing repeated locus from xtx.bf")
     xtx.bf <- xtx.bf[-which.in.xtx,] 
     REMOVE ALSO FROM locinames
    }
  }

# convert BF to log-BF, average rhos
	bfcols <- grep(".BF", colnames(xtx.bf))
  rho1cols <- grep(".rho1", colnames(xtx.bf))
  rho2cols <- grep(".rho2", colnames(xtx.bf))
   
  BFlog <- log(xtx.bf[,bfcols])
  rhoAve <- (xtx.bf[,rho1cols] + xtx.bf[,rho2cols])/2
  xtx.bf2 <- data.frame(xtx.bf[,1:3], BFlog, rhoAve)
  #head(xtx.bf2)

# reorder the matrices
  matches <- match(loci.names$Name, xtx.bf2$Name1)
  #xtx.bf2 <- matrix(NA, nrow=nrow(loci.names), ncol=ncol(xtx.bf)) 
  #xtx.bf2 <- xtx.bf[matches,]	
# combine matrices
  MAF.df <- read.table(MAF, header=TRUE)
  myout <- data.frame(loci.names, MAF.df, xtx.bf2[matches, ])
  print(head(myout[,1:10]))
  print(tail(myout[,1:10]))
  
# check they lined up right
  "Is the loci names lined up with XTX?"
  print(identical(myout$Name[!is.na(myout$Name1)], myout$Name1[!is.na(myout$Name1)]))  

# remove extras  
  rm(xtx.bf, xtx.bf2)
  rmcols <- c(which(colnames(myout)=="Name1"), which(colnames(myout)=="Name2"))

# add name column for Sam
  splitname <- strsplit(myout$Name,split = "__")
    #not sure if this will work...
  if (SEQCAP==TRUE){
   splitname2 <- matrix(unlist(splitname), ncol=3, byrow=TRUE)
     colnames(splitname2) <- c("TranscriptCont", "GenomicCont", "bp")
     NameContig <- paste(splitname2[,2], splitname2[,3], sep="__")
  }
  if (GBS==TRUE){
    splitname2 <- matrix(unlist(splitname), ncol=2, byrow=TRUE)
      colnames(splitname2) <- c("GenomicCont", "bp")
      NameContig <- myout$Name
  }

  myout2 <- data.frame(NameContig, splitname2, myout[,-rmcols])
  rownames(myout2) <- NULL
  head(myout2[,1:10])
  tail(myout2[,1:10])
# write myout
  write.table(myout2, paste(XTX.BF,"MAF.v2", sep=""), col.names=TRUE, row.names=FALSE)
  
# get cutoffs function
  get.cutoffs <- function(var){
    sort.BF<- sort(var)
    #hist(sort.BF)
    n <- length(sort.BF)
    max.V <- sort.BF[n]
    V.999 <- sort.BF[n*.999]
    V.99 <- sort.BF[n*.99]
    V.95 <- sort.BF[n*.95]
    V.90 <- sort.BF[n*.90]
    return(data.frame(V.999, V.99, V.95, V.90, max.V))
  }
  
# write recommended BF cutoff
  BFcols <- grep(".BF", colnames(myout2))
  BF.cutoffs <- get.cutoffs(as.vector(as.matrix(myout2[,BFcols])))
  colnames(BF.cutoffs) <- sub("V", "BF", colnames(BF.cutoffs))
  print(BF.cutoffs)
  
# write recommended XTX cutoff
  XTXcol <- grep("XTX", colnames(myout2))
  XTX.cutoffs <- get.cutoffs(myout2[,XTXcol])
  colnames(XTX.cutoffs) <- sub("V", "XTX", colnames(XTX.cutoffs))
  print(XTX.cutoffs)
  
 write.table(BF.cutoffs, paste(XTX.BF,"MAF.v2.BFcutoffs", sep=""), col.names=TRUE, row.names=FALSE)
  write.table(XTX.cutoffs, paste(XTX.BF,"MAF.v2.XTXcutoffs", sep=""), col.names=TRUE, row.names=FALSE)
}  