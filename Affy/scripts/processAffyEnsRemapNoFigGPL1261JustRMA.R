setwd("../Data")

library(affy)
library(affydata)
library(vsn)
library(mouse4302mmensgcdf) #### hgu133plus2hsensgcdf



Data <- ReadAffy(cdfname="mouse4302mmensgcdf") #####

normalize.AffyBatch.methods <- c(normalize.AffyBatch.methods,"vsn")

#eset <- expresso(Data,normalize.method="qspline",bgcorrect.method="rma",pmcorrect.method="pmonly",summary.method="liwong")


es2 = justRMA(cdfname="mouse4302mmensgcdf") ####

setwd("..")
dir.create("processed")
setwd("processed")


write.exprs(es2,file="rma.txt")
print("rma written")




save.image("Rdata")