# Written by Tyler William H Backman in 2016

# Algorithm BicBin
bicBinCluster <- function(binaryMatrix, totalClusters, alpha=0.5, beta=0.5, minDensityP1=0, tries=900){
    clusters <- list()
    for(i in 1:totalClusters){
        if(sum(binaryMatrix) == 0)
            break
        bicluster <- .algorithmA_2(binaryMatrix, alpha, beta, minDensityP1, tries)
        binaryMatrix[bicluster$rows, bicluster$cols] <- 0 # zero out this clusters region for next round
        clusters[[i]] <- list(rows=bicluster$rows, cols=bicluster$cols, score=bicluster$score)
    }
    return(clusters)
}

# Algorithm A_2 find clusters above set density
.algorithmA_2 <- function(binaryMatrix, alpha, beta, minDensityP1, tries){
    density <- -1 
    rowCoords <- 1:nrow(binaryMatrix)
    colCoords <- 1:ncol(binaryMatrix)
    while(density < minDensityP1){
        if(density < 0)
            density <- sum((binaryMatrix[rowCoords,colCoords])/(length(rowCoords) * length(colCoords))) 
        resultCoords <- .algorithmA_1(binaryMatrix[rowCoords, colCoords], alpha, beta, tries, density)
        rowCoords <- rowCoords[resultCoords$rows]
        colCoords <- colCoords[resultCoords$cols]
        density <- sum((binaryMatrix[rowCoords,colCoords])/(length(rowCoords) * length(colCoords))) 
    }
    return(list(rows=rowCoords, cols=colCoords, score=resultCoords$score))
}

# Algorithm A_1 in parallel
.algorithmA_1 <- function(binaryMatrix, alpha, beta, tries, p){
    allResults <- foreach(thisTry = 1:tries) %dorng% {
        res1 <- .BicBin(binaryMatrix, alpha, beta, p, proc_genes=TRUE) # Algorithm A_1 by columns
         res2 <- .BicBin(binaryMatrix, alpha, beta, p, proc_genes=FALSE) # Algorithm A_1 by rows
         if(res1$score > res2$score)
             return(res1)
         else
            return(res2)
    }
    maxCoords <- allResults[[which.max(sapply(allResults, function(x) x$score))]]
    return(list(rows=(maxCoords$x == 1), cols=(maxCoords$y == 1), score=maxCoords$score))
}
