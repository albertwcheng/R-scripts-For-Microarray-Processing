# Process Agilent text files for one-color-array and do between array normalization.
# produce G-norm files
# Albert Cheng
#
#
# Make sure the only ".txt" files in the directory are Agilent data files (or the script will crash)
#
# Put this is the same directory as the Agilent data files and run the command
#
# R --vanilla < normalize_summarize_Agilent_arrays_OneColor_withTargetFiles.R
#

origPath=getwd()
# Go to directory with *.txt files from scanner
setwd("../raw")


spotTypeFileName="SpotTypes.txt"
targetFileName="targets.txt"
outputPrefix="../tlimmaOutput"
withControlPath="withControl"

# This package needed
library(limma)

targets=readTargets(targetFileName)
spottypes <- readSpotTypes(spotTypeFileName)
#scanFiles = dir(pattern = ".*.txt$") #not using simple txt pattern match, use target file!~
num.arrays = length(targets$FileName)

# Only reads in the G values, fake the R values as G as well

maData = read.maimages(targets$FileName,names=targets$Cy3,
             columns = list(G = "gMeanSignal", Gb = "gBGUsed", R =
 "gProcessedSignal",
              Rb = "gBGMedianSignal"),
             annotation= c("Row", "Col", "FeatureNum", "ProbeUID",
 "ControlType",
              "ProbeName", "GeneName", "SystematicName"))


print("done reading")
maData$genes$Status <- controlStatus(spottypes,maData)


# Go back to starting directory
unlink(outputPrefix,recursive=TRUE)


dir.create(outputPrefix)
# Go to output directory
setwd(outputPrefix)
dir.create(withControlPath)

setwd(withControlPath)

# Do not background subtract (Default Agilent method ???) but add offset of 50
maData.nobg.0 = backgroundCorrect(maData, method="none", offset=0)
maData.nobg.50 = backgroundCorrect(maData, method="none", offset=50)


#G values only. Use Quantile Normalization for G channel
MAq0 = normalizeBetweenArrays(maData.nobg.0, method="Gquantile")
MAq50 = normalizeBetweenArrays(maData.nobg.50, method="Gquantile")



# Optional QC
# Densities of arrays before normalization
pdf("plotDensities_raw.pdf", h=8.5, w=11)
plotDensities(maData, singlechannels=( (ncol(maData$R)+1):(ncol(maData$R)+ncol(maData$G)) )) #plot only the green channel
dev.off()

# Densities of arrays after background correction
pdf("plotDensities_bgcorrected.0.pdf", h=8.5, w=11)
plotDensities(maData.nobg.0, singlechannels=( (ncol(maData.nobg.0$R)+1):(ncol(maData.nobg.0$R)+ncol(maData.nobg.0$G)) )) #plot only the green channel
dev.off()

pdf("plotDensities_bgcorrected.50.pdf", h=8.5, w=11)
plotDensities(maData.nobg.50, singlechannels=( (ncol(maData.nobg.50$R)+1):(ncol(maData.nobg.50$R)+ncol(maData.nobg.50$G)) )) #plot only the green channel
dev.off()

normRGq0=RG.MA(MAq0)
normRGq50=RG.MA(MAq50)


# Densities of arrays after quantile normalization between arrays
pdf("plotDensities_after_between_norm.bgc0.pdf", h=8.5, w=11)
plotDensities(normRGq0,  singlechannels=( (ncol(maData$R)+1):(ncol(maData$R)+ncol(maData$G)) ))
plotDensities(normRGq50,  singlechannels=( (ncol(maData$R)+1):(ncol(maData$R)+ncol(maData$G)) ))
dev.off()

# Print MA plots: raw; after 1st norm step; after 2nd norm step

#for (i in 1:num.arrays)
#{
#	png.file = paste("MA_plots_maData_", i, ".png", sep="")
#	png(png.file, h=1200, w=1600)
#	par(mfrow=c(2,3), mai=c(0.5,0.5,1,0.5) )		# c(bottom, left, top, right)
	
#	plotMA(maData, array=i, main=paste(colnames(maData)[i], " (raw)"))
#	plotMA(maData.nobg.50, array=i, main=" raw; offset=50")
#	plotMA(MA.loess.0, array=i, main="loess; offset=0")
#	plotMA(MA.loess.50, array=i, main="loess; offset=50")
#	plotMA(MA.loess.q.0, array=i, main="loess + Aquantile; offset=0")
#	plotMA(MA.loess.q.50, array=i, main="loess + Aquantile; offset=50")
#
#	dev.off()
#}




#write tables
write.table(cbind(maData$genes, log2(maData$G)), file="G.raw.txt", sep="\t", quote=FALSE, row.names=F)
write.table(cbind(maData.nobg.0$genes, log2(maData.nobg.0$G)), file="G.bgc0.txt", sep="\t", quote=FALSE, row.names=F)
write.table(cbind(normRGq0$genes, log2(normRGq0$G)), file="G.bgc0.norm.txt", sep="\t", quote=FALSE, row.names=F)
write.table(cbind(GeneName=normRGq0$genes$GeneName, log2(normRGq0$G)), file="G.bgc0.norm.gn.txt", sep="\t", quote=FALSE, row.names=F)

#go back to wherever we start off.
setwd(origPath)



################ EXIT . NO EXECUTION BEYOND THIS POINT ##########################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
