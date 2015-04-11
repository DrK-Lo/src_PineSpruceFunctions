### Assumes that some code from ProcessBFXTX_Spruce.R is in memory
#cd /data/seqcap/spruce/bwa_pseudo/final_tables/bayenv2Results10nonimputed_rm342479

MultiD2 <- function(df.vars, outfile, locusname){
# df.vars is a dataframe with each row a locus and each column a statsitic
# Calculate Malanahobis distance
  mu <- t(t(colMeans(df.vars, na.rm=TRUE))) # make it a column vector
  S <- cov(df.vars, use="pairwise.complete.obs") #variance-covariance matrix
  output <- 100000
  
  D <- NULL
	for (i in 1:nrow(df.vars)){
    x <- t(df.vars[i,])
    diff <- x-mu
    D[i] <- sqrt( t(diff) %*% solve(S) %*% (diff) )
    if(i==output){
      write.table(data.frame(D=D[1:output], Name2=locusname[1:output]), file=outfile, 
                  append=FALSE,row.names=FALSE, col.names=TRUE)
    }
    if(i>output & i%%output==0){
      write.table(data.frame(D[(i-output+1):i], locusname[(i-output+1):i]), file=outfile, 
                  append=TRUE,row.names=FALSE, col.names=FALSE)
    }
    if (i==nrow(df.vars)){
      last <- i%%output
      write.table(data.frame(D[(i-last+1):i], locusname[(i-last+1):i]), file=outfile, 
                  append=TRUE,row.names=FALSE, col.names=FALSE)
    }
  }
  return(D)
}


runGroupSpruce <- function(groupNum, file){ 
	setwd("/data/seqcap/spruce/bwa_pseudo/final_tables/bayenv2Results10nonimputed_rm342479")
	df <- read.table(file, header=TRUE)
	groups <- read.table("seqcap.bayenv.BF.envi2.clusters", header=TRUE)
	thisgroup <- paste(groups$enviAbb[groups$G6==groupNum], ".BF", sep="")
	column <- rep(NA, length(thisgroup))
	for (i in 1:length(thisgroup)){
		column[i] <- which(colnames(df)==thisgroup[i])
	}
	df.sub <- df[,column]
	outfile <- paste(file, ".D.G", groupNum, sep="")
	MultiD2(df.sub, outfile, df$Name)
}

runGroupPine <- function(groupNum, file){ 
	setwd("/data/seqcap/pine/bwa_pseudo/final_tables/bayenv2Results10nonimputed_rm580plus")
	df <- read.table(file, header=TRUE)
	groups <- read.table("var_out_pine_all_COMBINED.table.contig_flt10.bayenv.envi2.clusters", header=TRUE)
	thisgroup <- paste(groups$enviAbb[groups$G6==groupNum], ".BF", sep="")
	column <- rep(NA, length(thisgroup))
	for (i in 1:length(thisgroup)){
		column[i] <- which(colnames(df)==thisgroup[i])
	}
	df.sub <- df[,column]
	outfile <- paste(file, ".D.G", groupNum, sep="")
	MultiD2(df.sub, outfile, df$Name)
}