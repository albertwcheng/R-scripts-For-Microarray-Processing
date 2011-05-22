# Process Agilent text files and do loess and Quantile (over all arrays and channels) normalization
#
# Modified by Albert Cheng - MIT CSBi
# awcheng@mit.edu
# Mar 28, 2009
#
#
#
# For more info about commands, see
# See http://www.bioconductor.org/packages/1.9/bioc/vignettes/limma/inst/doc/usersguide.pdf
# which is a User's Guide for the "limma" package in R
#
# Make sure the only ".txt" files in the directory are Agilent data files (or the script will crash)
#
# Put this is the same directory as the Agilent data files and run the command
#
# /nfs/BaRC/R/R_vanilla < normalize_summarize_Agilent_arrays.R
#
# or to run it on the LSF cluster,
#
# bsub "/nfs/BaRC/R/R_vanilla < normalize_summarize_Agilent_arrays.R"
#

# Go to directory with *.txt files from scanner

source("arrayTable.R")

outputPrefix="../tlimmaOutput"

setwd(outputPrefix)

load("ns.RData")

setwd(withControlPath)

# This package needed
library(limma)

MA.loess.q.0 = normalizeBetweenArrays(MA.loess.0, method="quantile")
MA.loess.q.50 = normalizeBetweenArrays(MA.loess.50, method="quantile")


# Densities of arrays after quantile normalization between arrays
pdf("plotDensities_after_loess_quantile_norm.pdf", h=8.5, w=11)
plotDensities(MA.loess.q.0)
plotDensities(MA.loess.q.50)
dev.off()

# Print MA plots: raw; after 1st norm step; after 2nd norm step

for (i in 1:num.arrays)
{
	png.file = paste("MA_plots_maData_quantile", i, ".png", sep="")
	png(png.file, h=1200, w=1600)
	par(mfrow=c(2,3), mai=c(0.5,0.5,1,0.5) )		# c(bottom, left, top, right)
	
	plotMA(maData, array=i, main=paste(colnames(maData)[i], " (raw)"))
	plotMA(maData.nobg.50, array=i, main=" raw; offset=50")
	plotMA(MA.loess.0, array=i, main="loess; offset=0")
	plotMA(MA.loess.50, array=i, main="loess; offset=50")
	plotMA(MA.loess.q.0, array=i, main="loess + quantile; offset=0")
	plotMA(MA.loess.q.50, array=i, main="loess + quantile; offset=50")

	dev.off()
}

# Print M and A values

# Before quantile normalization
MA.loess.0.NT=getNameTransformedMA(MA.loess.0,targets)
writeMATable(MA.loess.0.NT, "log2Ratio_loess.txt")
# After quantile normalization
MA.loess.q.0.NT=getNameTransformedMA(MA.loess.q.0,targets)
writeMATable(MA.loess.q.0.NT,"log2Ratio_loess_q.txt")


png("boxplot_raw.png", h=1200, w=1600)
par(las=3, omi=c(0.5, 0, 0, 0), cex=0.75) 
boxplot( as.matrix(maData.nobg.0) ~col(as.matrix(maData.nobg.0)), notch=T, main = "Rawl; no BG correction")
dev.off()

png("boxplot_loess.png", h=1200, w=1600)
par(las=3, omi=c(0.5, 0, 0, 0), cex=0.75) 
boxplot( as.matrix(MA.loess.0) ~col(as.matrix(MA.loess.0)), notch=T, main = "With loess normalization")
dev.off()

png("boxplot_loess_q.png", h=1200, w=1600)
par(las=3, omi=c(0.5, 0, 0, 0), cex=0.75) 
boxplot( as.matrix(MA.loess.q.0) ~col(as.matrix(MA.loess.q.0)), notch=T, main = "With loess + Aquantile normalization")
dev.off()


# Convert MA data into RG data
RG.loess.0=RG.MA(MA.loess.0)
RG.loess.0.NT.log2=log2RG(getNameTransformedRG(RG.loess.0, targets))
writeRGTable(RG.loess.0.NT.log2,"log2Intensity_loess.txt")

RG.loess.q.0  = RG.MA(MA.loess.q.0)
RG.loess.q.0.NT.log2= log2RG(getNameTransformedRG(RG.loess.q.0,targets))
writeRGTable(RG.loess.q.0.NT.log2,"log2Intensity_loess_q.txt")

#now write everything!
writeRGMATable(RG.loess.0.NT.log2,MA.loess.0.NT,"expSummary_loess.txt")
writeRGMATable(RG.loess.q.0.NT.log2,MA.loess.q.0.NT,"expSummary_loess_q.txt")

setwd("..")

saveImageWithBackup()



