setwd("../Data")

library(affy)
library(affydata)
library(vsn)
library(mogene10stv1mmensgcdf) #### hgu133plus2hsensgcdf


Data <- ReadAffy(cdfname="mogene10stv1mmensgcdf") #####

es2 = justRMA(cdfname="mouse4302mmensgcdf") ####

setwd("..")
dir.create("processed")
setwd("processed")


write.exprs(es2,file="rma.mogene10stv1mmensgcdf.txt")
print("rma written")




save.image("Rdata")