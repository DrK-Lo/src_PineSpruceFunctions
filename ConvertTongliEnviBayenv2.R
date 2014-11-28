setwd("/Users/katie/Desktop/AdaptreeData")

sp.envi <- read.csv("var_out_588_env.csv")
head(sp.envi)

which(colnames(sp.envi)=="Working.SL")

remove.cols <- c(1,2,3,4,27)

rel.envi <- sp.envi[,-remove.cols] #relevant columns

std.envi <- rel.envi #will replace with standardized values
tail(std.envi)

for (i in 1:ncol(rel.envi)){
  colm <- rel.envi[,i]
  colm[which(colm==-9999)]=NA
  print(i)
    if(sum(is.na(colm)>0)){print(which(is.na(colm)))}
  remove <- which(is.na(colm))
  std.envi[,i]<-colm-mean(colm, na.rm=TRUE)
}

head(std.envi)
tail(std.envi)
write.table(x=colnames(std.envi), "EnvironmentNames.txt", ncol=1, col.names=FALSE, row.names=FALSE)