setwd("../Data")

library(affy)
library(affydata)
library(vsn)
library(mouse4302mmensgcdf) #### hgu133plus2hsensgcdf


Data <- ReadAffy(cdfname="mouse4302mmensgcdf") #####

es2 = justRMA(cdfname="mouse4302mmensgcdf") ####

setwd("..")
dir.create("processed")
setwd("processed")


write.exprs(es2,file="rma.mouse4302mmensgcdf.txt")
print("rma written")




save.image("Rdata.mouse4302mmensgcdf.justRMA")