setwd("../Data")

library(affy)
library(affydata)
library(vsn)
library(hgu133plus2hsensgcdf) #### hgu133plus2hsensgcdf



Data <- ReadAffy(cdfname="hgu133plus2hsensgcdf") #####

es2 = justRMA(cdfname="hgu133plus2hsensgcdf") ####

setwd("..")
dir.create("processed")
setwd("processed")


write.exprs(es2,file="rma.hgu133plus2hsensgcdf.txt")
print("rma written")

es3 = mas5(Data)

write.exprs(es3,file="mas5.hgu133plus2hsensgcdf.txt")
print("mas5 written")

esMasCall = mas5calls(Data)

write.exprs(esMasCall,file="mas5calls.hgu133plus2hsensgcdf.txt")
print("mas5calls written")

es1 = expresso(Data,bg.correct=FALSE, normalize.method="vsn",pmcorrect.method="pmonly",summary.method="medianpolish")

write.exprs(es1,file="vsn.hgu133plus2hsensgcdf.txt")


save.image("Rdata.hgu133plus2hsensgcdf")