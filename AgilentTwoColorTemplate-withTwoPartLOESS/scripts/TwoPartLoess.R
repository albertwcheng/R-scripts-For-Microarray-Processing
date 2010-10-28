library(limma)

doLOESS <- function(maData,threshold)
{
	if(threshold<0)
	{
		return(normalizeWithinArrays(maData, method="loess"))
	}
	else
	{
		return(twoPartLoess(maData,threshold))
	}
}

twoPartLoess <- function(maData,threshold)
{
	#first do a overall LOESS
	
	MA.loess=normalizeWithinArrays(maData,method="loess",bc.method="none")
	
	num.arrays=dim(MA.loess$A)[2]
	MA.twoPartLoess=MA.loess #copy
	
	for( i in 1:num.arrays){
		this.array=MA.loess[,i]
		high.int=this.array$A>=threshold
		other.int=this.array$A<threshold
		this.array.high.int=this.array[high.int,]
		this.array.other.int=this.array[other.int,]
		this.array.high.int.loess=normalizeWithinArrays(this.array.high.int,method="loess",bc.method="none")
		this.array.other.int.loess=normalizeWithinArrays(this.array.other.int,method="loess",bc.method="none") #non changed to none June2, 2009
		MA.loess.merged=rbind(this.array.other.int.loess,this.array.high.int.loess)
		
		MA.loess.merged.sorted=MA.loess.merged[order(MA.loess.merged$genes[,1],MA.loess.merged$genes[,2]),]
		
		MA.twoPartLoess$M[,i]=MA.loess.merged.sorted$M
		MA.twoPartLoess$A[,i]=MA.loess.merged.sorted$A
		MA.twoPartLoess$genes=MA.loess.merged.sorted$genes
		
		if(!identical(MA.loess$genes[1:100000,],MA.twoPartLoess$genes[1:100000,]))
		{
			cat("error: After twoPartLoess() Gene Info Not Aligned")
			quit("no")
		}
	
	}
	

	
	return (MA.twoPartLoess)
	
}