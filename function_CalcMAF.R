#get MAF cond
# spruce:  ssh klott@rogue.zoology.ubc.ca:/data/seqcap/spruce/bwa_pseudo/final_tables 
  # screen; R; MakeMAFdf("var_out_588ind_all_v2_COMBINED.table.contig_flt10.bayenv")
  # screen; R; MakeMAFdf("spruce_GBS_all.table.contig_flt10.bayenv")
# pine: ssh klott@hulk.zoology.ubc.ca:/data/seqcap/pine/bwa_pseudo/final_tables
  # screen; R; MakeMAFdf("var_out_pine_all_COMBINED.table.contig_flt10.bayenv")
  # screen; R; MakeMAFdf("pine_GBS_all.table.contig_flt10.bayenv")

MakeMAFdf <- function(filename){
  be.df <- read.table(filename)
  nlines <- nrow(be.df)
  even <- seq(2, nlines, 2)
  odd <- seq(1, nlines,2)

  A1count <- rowSums(be.df[odd,])
  A2count <- rowSums(be.df[even,])

  A1freq <- A1count/(A1count + A2count)

  MAF <- A1freq
  cond <- A1freq > 0.5
  MAF[cond] <- 1 - A1freq[cond]

#hist(MAF)
  MAF <- as.numeric(MAF)
  He <- as.numeric(2*A1freq*(1-A1freq))
  He.gt.05 <- (He >= 0.05)
  print(sum(He.gt.05))
  write.table(data.frame(MAF, He, He.gt.05), file = paste(filename,"MAF", sep="."), col.names=TRUE, row.names=FALSE)
}