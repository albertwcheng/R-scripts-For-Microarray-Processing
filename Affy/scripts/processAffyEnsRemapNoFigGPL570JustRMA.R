setwd("../Data")

library(affy)
library(affydata)
library(vsn)
library(hgu133plus2hsensgcdf) #### hgu133plus2hsensgcdf



Data <- ReadAffy(cdfname="hgu133plus2hsensgcdf") #####

normalize.AffyBatch.methods <- c(normalize.AffyBatch.methods,"vsn")

#eset <- expresso(Data,normalize.method="qspline",bgcorrect.method="rma",pmcorrect.method="pmonly",summary.method="liwong")


es2 = justRMA(cdfname="hgu133plus2hsensgcdf") ####

setwd("..")
dir.create("processed")
setwd("processed")


write.exprs(es2,file="rma.txt")
print("rma written")

save.image("Rdata")