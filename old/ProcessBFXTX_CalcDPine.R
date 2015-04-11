
### Seqcap Pine: calculate multivariate BF for each of 6 groups
>screen
>R
source("ProcessBFXTX_CalcD.R")  
runGroupPine(1, "Pine_SeqCap.XTX.BFMAF.v2")

>screen
>R
source("ProcessBFXTX_CalcD.R")  
runGroupPine(2, "Pine_SeqCap.XTX.BFMAF.v2")

>screen
>R
source("ProcessBFXTX_CalcD.R")  
runGroupPine(3, "Pine_SeqCap.XTX.BFMAF.v2")


>screen
>R
source("ProcessBFXTX_CalcD.R")  
runGroupPine(4, "Pine_SeqCap.XTX.BFMAF.v2")


>screen
>R
source("ProcessBFXTX_CalcD.R")  
runGroupPine(5, "Pine_SeqCap.XTX.BFMAF.v2")

>screen
>R
source("ProcessBFXTX_CalcD.R")  
runGroupPine(6, "Pine_SeqCap.XTX.BFMAF.v2")


### Seqcap GBS: calculate multivariate BF for each of 6 groups
>screen
>R
file <- "Pine_GBS.XTX.BFMAF.v2"
  source("ProcessBFXTX_CalcD.R")  
runGroupPine(file, 1)

>screen
>R
file <- "Pine_GBS.XTX.BFMAF.v2"
  source("ProcessBFXTX_CalcD.R")  
runGroupPine(file, 2)

>screen
>R
file <- "Pine_GBS.XTX.BFMAF.v2"
  source("ProcessBFXTX_CalcD.R")  
runGroupPine(file, 3)

>screen
>R
file <- "Pine_GBS.XTX.BFMAF.v2"
  source("ProcessBFXTX_CalcD.R")  
runGroupPine(file, 4)

>screen
>R
file <- "Pine_GBS.XTX.BFMAF.v2"
  source("ProcessBFXTX_CalcD.R")  
runGroupPine(file, 5)

>screen
>R
file <- "Pine_GBS.XTX.BFMAF.v2"
  source("ProcessBFXTX_CalcD.R")  
runGroupPine(file, 6)
