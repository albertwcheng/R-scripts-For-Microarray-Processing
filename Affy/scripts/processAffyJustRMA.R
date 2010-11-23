setwd("../Data")

library(affy)
library(affydata)
library(vsn)



Data <- ReadAffy() #####

es2 = justRMA() ####

setwd("..")
dir.create("processed")
setwd("processed")


write.exprs(es2,file="rma.txt")
print("rma written")

save.image("Rdata")