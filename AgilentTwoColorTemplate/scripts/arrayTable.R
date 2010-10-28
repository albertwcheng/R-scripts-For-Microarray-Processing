# Swap in proper sample names for RG and MA list
# Albert Cheng - MIT CSBi
# awcheng@mit.edu
# Mar 28, 2009
#
# Required columns for targets
# Cy5  -  name of the Cy5 sample
# Cy3  -  name of the Cy3 sample
# M	   -  positive means ratio should be Cy5-Cy3, negative means ratio be Cy3-Cy5



"%dm%"<-function(X,Y){
	for(i in 1:length(Y)){X[,i]<-X[,i]*Y[i]} 
	return(X)
	
	}

swapMatrixColumns<-function(X,Y,swapMatrix){
	nX<-X
	nY<-Y
	for(i in 1:length(swapMatrix)){
			if( swapMatrix[i] < 0){
				nX[,i]<-Y[,i]
				nY[,i]<-X[,i]
			}
		}
	
	return(list(X=nX,Y=nY))
	
	}

swapArray <-function(X,Y,swapMatrix){
	nX<-X
	nY<-Y
	for(i in 1:length(swapMatrix)){
			if( swapMatrix[i] < 0){
				nX[i]<-Y[i]
				nY[i]<-X[i]
			}
		}
	
	return(list(X=nX,Y=nY))	
	
	}

getNameTransformedRG<-function(RG,targets){
	newRG<-RG
	L<-swapMatrixColumns(RG$R,RG$G,targets$M)
	newRG$R<-L$X
	newRG$G<-L$Y
	L<-swapArray(targets$Cy5,targets$Cy3,targets$M)
	Rnames<-L$X
	Gnames<-L$Y
	colnames(newRG$R)<-Rnames
	colnames(newRG$G)<-Gnames
	
	return(newRG)
	}
	
getNameTransformedMA<-function(MA,targets){
	newMA<-MA
	newMA$M<-newMA$M %dm% targets$M
	
	L<-swapArray(targets$Cy5,targets$Cy3,targets$M)
	Rnames<-L$X
	Gnames<-L$Y
	Mnames<-paste(Rnames,Gnames,sep="-")
	Anames<-paste("mean(",paste(Rnames,Gnames,sep=","),")",sep="")
	colnames(newMA$M)<-Mnames
	colnames(newMA$A)<-Anames
	return(newMA)
	
	}

writeMTable<-function(MA,filename){
	write.table(cbind(MA$genes, MA$M), file=filename, sep="\t", quote=FALSE, row.names=F)
	}
	
writeATable<-function(MA,filename){
	write.table(cbind(MA$genes, MA$A), file=filename, sep="\t", quote=FALSE, row.names=F)
	}
	
writeMATable<-function(MA,filename){
	write.table(cbind(MA$genes, MA$A,MA$M), file=filename, sep="\t", quote=FALSE, row.names=F)

	}
	
writeRTable<-function(RG,filename){
	write.table(cbind(RG$genes, RG$R), file=filename, sep="\t", quote=FALSE, row.names=F)
	}
	
writeGTable<-function(RG,filename){
	write.table(cbind(RG$genes, RG$G), file=filename, sep="\t", quote=FALSE, row.names=F)
	}
	
writeRGTable<-function(RG,filename){
	write.table(cbind(RG$genes,RG$G,RG$R), file=filename, sep="\t", quote=FALSE, row.names=F)
	}
	
writeRGTable2<-function(RG, filename ,extra){


		write.table(cbind(RG$genes,RG$G,RG$R,extra), file=filename, sep="\t", quote=FALSE, row.names=F)
		

}
		
writeRGMATable<-function(RG,MA,filename){
	write.table(cbind(RG$genes,RG$G,RG$R,MA$A,MA$M),file=filename,sep="\t",quote=FALSE,row.names=F)
	}
	
writeRGMATable2 <-function(RG,MA, extra,filename){
	write.table(cbind(RG$genes,RG$G,RG$R,MA$A,MA$M, extra),file=filename,sep="\t",quote=FALSE,row.names=F)
	}
	
log2RG<-function(RG){
	RG$R=log2(RG$R)
	RG$G=log2(RG$G)
	return (RG)
	}
	
saveImageWithBackup<-function(rmTmpV=TRUE,pattern="_00000"){
	timestamp=gsub(" ","_",gsub(":"," ",date()))
	if(rmTmpV){
			rmTmpVars(pattern)
		}
	file.rename("ns.RData",paste("ns.RData.bk_",timestamp))
	save.image(file="ns.RData")
	
	
	}
	
rmTmpVars <- function(patternOfVar="_00000"){

		varMatched <- ls(globalenv(),pattern=patternOfVar)
		print("variables to kill:")
		print(varMatched)
		rm(list=varMatched, envir= globalenv())

	}