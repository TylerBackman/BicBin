# (C) 2016 Tyler William H Backman

# perform bicbin clustering in parallel
bicBinCluster <- function(binaryMatrix, totalClusters, cores, alpha=0.5, beta=0.5, minDensityP1=0, tries=900){
    M <- nrow(binaryMatrix)
    N <- ncol(binaryMatrix)
    p <- sum(binaryMatrix) / (M*N)
    clusters <- list()
    for(i in 1:totalClusters){
	density <- 0
	while(density < minDensityP1){
        	allResults <- foreach(thisTry = 1:tries) %dorng% {
             	 res1 <- .BicBin(binaryMatrix, alpha, beta, p, proc_genes=TRUE) # by columns
             	 res2 <- .BicBin(binaryMatrix, alpha, beta, p, proc_genes=FALSE) # by rows
            	  if(res1$score > res2$score)
            	    return(res1)
            	  else
             	   return(res2)
        	}
		max <- allResults[[which.max(sapply(allResults, function(x) x$score))]]
		subMatrix <- binaryMatrix[max$x == 1, max$y == 1]
		density <- sum(subMatrix) / (nrow(subMatrix) * ncol(subMatrix))
	}
        binaryMatrix[max$x == 1, max$y == 1] <- 0 # zero out this clusters region for next round
        clusters[[i]] <- list(rows=row.names(binaryMatrix)[max$x == 1], 
		cols=colnames(binaryMatrix)[max$y == 1], score=max$score)
    }
    return(clusters)
}
