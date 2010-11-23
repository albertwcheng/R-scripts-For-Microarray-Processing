setwd("../Data")

library(affy)
library(affydata)
library(vsn)
library(mouse4302mmensgcdf) #### hgu133plus2hsensgcdf



Data <- ReadAffy(cdfname="mouse4302mmensgcdf") #####

normalize.AffyBatch.methods <- c(normalize.AffyBatch.methods,"vsn")


es2 = justRMA(cdfname="mouse4302mmensgcdf") ####

setwd("..")
dir.create("processed")
setwd("processed")


write.exprs(es2,file="rma.mouse4302mmensgcdf.txt")
print("rma written")

es3 = mas5(Data)

write.exprs(es3,file="mas5.mouse4302mmensgcdf.txt")
print("mas5 written")

esMasCall = mas5calls(Data)

write.exprs(esMasCall,file="mas5calls.mouse4302mmensgcdf.txt")
print("mas5calls written")

es1 = expresso(Data,bg.correct=FALSE, normalize.method="vsn",pmcorrect.method="pmonly",summary.method="medianpolish")

write.exprs(es1,file="vsn.mouse4302mmensgcdf.txt")


save.image("Rdata")