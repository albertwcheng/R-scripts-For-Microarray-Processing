setwd("../Data")

library(affy)
library(affydata)
library(vsn)
library(mogene10stv1mmensgcdf) #### hgu133plus2hsensgcdf


Data <- ReadAffy(cdfname="mogene10stv1mmensgcdf") #####

#normalize.AffyBatch.methods <- c(normalize.AffyBatch.methods,"vsn")

es2 = justRMA(cdfname="mouse4302mmensgcdf") ####

setwd("..")
dir.create("processed")
setwd("processed")


write.exprs(es2,file="rma.txt")
print("rma written")




save.image("Rdata")