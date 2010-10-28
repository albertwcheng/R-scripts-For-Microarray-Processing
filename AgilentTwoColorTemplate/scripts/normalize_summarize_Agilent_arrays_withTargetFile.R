# Process Agilent text files and do loess and Aquantile normalization
#
# Modified by Albert Cheng - MIT CSBi
# awcheng@mit.edu
# Mar 28, 2009
#
#

source("arrayTable.R")

origPath=getwd()

setwd("../raw")





spotTypeFileName="SpotTypes.txt"
targetFileName="targets.txt"
outputPrefix="../tlimmaOutput"
withControlPath="withControl"
noControlPath="noControl"

# This package needed
library(limma)

targets=readTargets(targetFileName)
spottypes <- readSpotTypes(spotTypeFileName)
#scanFiles = dir(pattern = ".*.txt$")
num.arrays = length(targets$FileName)

# This program reads gMeanSignal and rMeanSignal (for foreground)
# and gBGMedianSignal and rBGMedianSignal (for background)
maData = read.maimages(targets, source="agilent")
maData$genes$Status <- controlStatus(spottypes,maData)


# read agilent file, not counting weight for control group
f = function(x) as.numeric(x$ControlType == 0)
maData.noControls = read.maimages(targets, source="agilent", wt.fun=f)
maData.noControls$genes$Status <- controlStatus(spottypes,maData.noControls)

#clean output dir
unlink(outputPrefix,recursive=TRUE)


dir.create(outputPrefix)
# Go to output directory
setwd(outputPrefix)
dir.create(withControlPath)
dir.create(noControlPath)


#save raw data as RData files #kinda redundant? remove!
#save(maData, file="raw_maData.Rdata")
#save(maData.noControls, file="raw_maData.noControls.Rdata")


setwd(withControlPath)


# Do not background subtract (Default Agilent method ???) but add offset of 50
maData.nobg.0 = backgroundCorrect(maData, method="none", offset=0)
maData.nobg.50 = backgroundCorrect(maData, method="none", offset=50)

# Normalize using loess (intensity-dependent)
MA.loess.0 = normalizeWithinArrays(maData.nobg.0, method="loess")
MA.loess.0$genes$Status <- controlStatus(spottypes,MA.loess.0)
MA.loess.50 = normalizeWithinArrays(maData.nobg.50, method="loess")
MA.loess.50$genes$Status <- controlStatus(spottypes,MA.loess.50)
# Normalize between arrays (using quantile method on average (A) intensities)
MA.loess.aq.0 = normalizeBetweenArrays(MA.loess.0, method="Aquantile")
MA.loess.aq.50 = normalizeBetweenArrays(MA.loess.50, method="Aquantile")


# Optional QC
# Densities of arrays before normalization
pdf("plotDensities_before_norm.pdf", h=8.5, w=11)
plotDensities(maData)
dev.off()
# Densities of arrays after loess normalization within arrays
pdf("plotDensities_after_loess_norm.pdf", h=8.5, w=11)
plotDensities(MA.loess.0)
plotDensities(MA.loess.50)
dev.off()
# Densities of arrays after quantile normalization between arrays
pdf("plotDensities_after_loess_aquantile_norm.pdf", h=8.5, w=11)
plotDensities(MA.loess.aq.0)
plotDensities(MA.loess.aq.50)
dev.off()

# Print MA plots: raw; after 1st norm step; after 2nd norm step

for (i in 1:num.arrays)
{
	png.file = paste("MA_plots_maData_aquantile_", i, ".png", sep="")
	png(png.file, h=1200, w=1600)
	par(mfrow=c(2,3), mai=c(0.5,0.5,1,0.5) )		# c(bottom, left, top, right)
	
	plotMA(maData, array=i, main=paste(colnames(maData)[i], " (raw)"))
	plotMA(maData.nobg.50, array=i, main=" raw; offset=50")
	plotMA(MA.loess.0, array=i, main="loess; offset=0")
	plotMA(MA.loess.50, array=i, main="loess; offset=50")
	plotMA(MA.loess.aq.0, array=i, main="loess + Aquantile; offset=0")
	plotMA(MA.loess.aq.50, array=i, main="loess + Aquantile; offset=50")

	dev.off()
}

# Print M and A values

# Before quantile normalization
MA.loess.0.NT=getNameTransformedMA(MA.loess.0,targets)
writeMATable(MA.loess.0.NT, "log2Ratio_loess.txt")
# After quantile normalization
MA.loess.aq.0.NT=getNameTransformedMA(MA.loess.aq.0,targets)
writeMATable(MA.loess.aq.0.NT,"log2Ratio_loess_q.txt")


####write.table(cbind(MA.loess.0$genes, MA.loess.0$M), file="M_log2ratios_loess.TXT", sep="\t", quote=FALSE, row.names=F)
####write.table(cbind(MA.loess.0$genes, MA.loess.0$A), file="A_log2intensities_loess.TXT", sep="\t", quote=FALSE, row.names=F)

#### After quantile normalization
####write.table(cbind(MA.loess.aq.0$genes, MA.loess.aq.0$M), file="M_log2ratios_loess_q.TXT", sep="\t", quote=FALSE, row.names=F)
####write.table(cbind(MA.loess.aq.0$genes, MA.loess.aq.0$A), file="A_log2intensities_loess_q.TXT", sep="\t", quote=FALSE, row.names=F)

png("boxplot_raw.png", h=1200, w=1600)
par(las=3, omi=c(0.5, 0, 0, 0), cex=0.75) 
boxplot( as.matrix(maData.nobg.0) ~col(as.matrix(maData.nobg.0)), notch=T, main = "Rawl; no BG correction")
dev.off()

png("boxplot_loess.png", h=1200, w=1600)
par(las=3, omi=c(0.5, 0, 0, 0), cex=0.75) 
boxplot( as.matrix(MA.loess.0) ~col(as.matrix(MA.loess.0)), notch=T, main = "With loess normalization")
dev.off()

png("boxplot_loess_aq.png", h=1200, w=1600)
par(las=3, omi=c(0.5, 0, 0, 0), cex=0.75) 
boxplot( as.matrix(MA.loess.aq.0) ~col(as.matrix(MA.loess.aq.0)), notch=T, main = "With loess + Aquantile normalization")
dev.off()


# Convert MA data into RG data
RG.loess.0=RG.MA(MA.loess.0)
RG.loess.0.NT.log2=log2RG(getNameTransformedRG(RG.loess.0, targets))
writeRGTable(RG.loess.0.NT.log2,"log2Intensity_loess.txt")

RG.loess.aq.0  = RG.MA(MA.loess.aq.0)
RG.loess.aq.0.NT.log2=log2RG(getNameTransformedRG(RG.loess.aq.0,targets))
writeRGTable(RG.loess.aq.0.NT.log2,"log2Intensity_loess_aq.txt")

#now write everything!
writeRGMATable(RG.loess.0.NT.log2,MA.loess.0.NT,"expSummary_loess.txt")
writeRGMATable(RG.loess.aq.0.NT.log2,MA.loess.aq.0.NT,"expSummary_loess_aq.txt")

##### Get R values (channel with experimental samples) after quantile normalization
#####write.table(cbind(RG.loess.aq.0$genes, RG.loess.aq.0$R), file="R_intensities_loess_q.TXT", sep="\t", quote=FALSE, row.names=F)



##############################################################
setwd("..")
setwd(noControlPath)


# Do not background subtract (Default Agilent method ???) but add offset of 50
maData.nobg.0.noControls = backgroundCorrect(maData.noControls, method="none", offset=0)
maData.nobg.50.noControls  = backgroundCorrect(maData.noControls, method="none", offset=50)

# Normalize using loess (intensity-dependent)
MA.loess.0.noControls  = normalizeWithinArrays(maData.nobg.0.noControls , method="loess")
MA.loess.0.noControls$genes$Status <- controlStatus(spottypes,MA.loess.0.noControls)
MA.loess.50.noControls  = normalizeWithinArrays(maData.nobg.50.noControls , method="loess")
MA.loess.50.noControls$genes$Status <- controlStatus(spottypes,MA.loess.50.noControls)
# Normalize between arrays (using quantile method on average (A) intensities)
MA.loess.aq.0.noControls  = normalizeBetweenArrays(MA.loess.0.noControls , method="Aquantile")
MA.loess.aq.50.noControls  = normalizeBetweenArrays(MA.loess.50.noControls , method="Aquantile")


# Optional QC
# Densities of arrays before normalization
pdf("plotDensities_before_norm.pdf", h=8.5, w=11)
plotDensities(maData.noControls)
dev.off()
# Densities of arrays after loess normalization within arrays
pdf("plotDensities_after_loess_norm.pdf", h=8.5, w=11)
plotDensities(MA.loess.0.noControls)
plotDensities(MA.loess.50.noControls)
dev.off()
# Densities of arrays after quantile normalization between arrays
pdf("plotDensities_after_loess_aquantile_norm.pdf", h=8.5, w=11)
plotDensities(MA.loess.aq.0.noControls)
plotDensities(MA.loess.aq.50.noControls)
dev.off()

# Print MA plots: raw; after 1st norm step; after 2nd norm step

for (i in 1:num.arrays)
{
	png.file = paste("MA_plots_maData_aquantile_", i, ".png", sep="")
	png(png.file, h=1200, w=1600)
	par(mfrow=c(2,3), mai=c(0.5,0.5,1,0.5) )		# c(bottom, left, top, right)
	
	plotMA(maData.noControls, array=i, main=paste(colnames(maData)[i], " (raw)"))
	plotMA(maData.nobg.50.noControls, array=i, main=" raw; offset=50")
	plotMA(MA.loess.0.noControls, array=i, main="loess; offset=0")
	plotMA(MA.loess.50.noControls, array=i, main="loess; offset=50")
	plotMA(MA.loess.aq.0.noControls, array=i, main="loess + Aquantile; offset=0")
	plotMA(MA.loess.aq.50.noControls, array=i, main="loess + Aquantile; offset=50")

	dev.off()
}

# Print M and A values

# Before quantile normalization
MA.loess.0.noControls.NT=getNameTransformedMA(MA.loess.0.noControls,targets)
writeMATable(MA.loess.0.noControls.NT, "log2Ratio_loess.txt")
# After quantile normalization
MA.loess.aq.0.noControls.NT=getNameTransformedMA(MA.loess.aq.0.noControls,targets)
writeMATable(MA.loess.aq.0.noControls.NT,"log2Ratio_loess_q.txt")


# Before quantile normalization
###write.table(cbind(MA.loess.0.noControls$genes, MA.loess.0.noControls$M), file="M_log2ratios_loess.TXT", sep="\t", quote=FALSE, row.names=F)
###write.table(cbind(MA.loess.0.noControls$genes, MA.loess.0.noControls$A), file="A_log2intensities_loess.TXT", sep="\t", quote=FALSE, row.names=F)

# After quantile normalization
###write.table(cbind(MA.loess.aq.0.noControls$genes, MA.loess.aq.0.noControls$M), file="M_log2ratios_loess_q.TXT", sep="\t", quote=FALSE, row.names=F)
###write.table(cbind(MA.loess.aq.0.noControls$genes, MA.loess.aq.0.noControls$A), file="A_log2intensities_loess_q.TXT", sep="\t", quote=FALSE, row.names=F)

png("boxplot_raw.png", h=1200, w=1600)
par(las=3, omi=c(0.5, 0, 0, 0), cex=0.75) 
boxplot( as.matrix(maData.nobg.0.noControls) ~col(as.matrix(maData.nobg.0.noControls)), notch=T, main = "Rawl; no BG correction")
dev.off()

png("boxplot_loess.png", h=1200, w=1600)
par(las=3, omi=c(0.5, 0, 0, 0), cex=0.75) 
boxplot( as.matrix(MA.loess.0.noControls) ~col(as.matrix(MA.loess.0.noControls)), notch=T, main = "With loess normalization")
dev.off()

png("boxplot_loess_aq.png", h=1200, w=1600)
par(las=3, omi=c(0.5, 0, 0, 0), cex=0.75) 
boxplot( as.matrix(MA.loess.aq.0.noControls) ~col(as.matrix(MA.loess.aq.0.noControls)), notch=T, main = "With loess + Aquantile normalization")
dev.off()


RG.loess.0.noControls=RG.MA(MA.loess.0.noControls)
RG.loess.0.noControls.NT.log2=log2RG(getNameTransformedRG(RG.loess.0.noControls,targets))
writeRGTable(RG.loess.0.noControls.NT.log2,"log2Intensity_loess.txt")

RG.loess.aq.0.noControls  = RG.MA(MA.loess.aq.0.noControls)
RG.loess.aq.0.noControls.NT.log2=log2RG(getNameTransformedRG(RG.loess.aq.0.noControls,targets))
writeRGTable(RG.loess.aq.0.noControls.NT.log2,"log2Intensity_loess_aq.txt")

#now write everything!
writeRGMATable(RG.loess.0.noControls.NT.log2,MA.loess.0.noControls.NT,"expSummary_loess.txt")
writeRGMATable(RG.loess.aq.0.noControls.NT.log2,MA.loess.aq.0.noControls.NT,"expSummary_loess_aq.txt")

setwd("..")

save.image(file="ns.RData")

setwd(origPath)