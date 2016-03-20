\name{bicBinCluster}
\alias{bicBinCluster}
\title{
	Biclustering sparse binary data in parallel
}
\description{
Biclusters sparse binary data in parallel across many cores based on BicBin by Miranda van Uitert, Wouter Meuleman, and Lodewyk Wessels. 
}
\usage{
bicBinCluster(binaryMatrix, totalClusters, alpha = 0.5, beta = 0.5, minDensityP1 = 0, tries = 900)
}
\arguments{
  \item{binaryMatrix}{
This is a binary matrix (0s and 1s) with unique row and column names.
}
  \item{totalClusters}{
The desired number of clusters to find.
}
  \item{alpha}{
The alpha term from the BicBin algorithm indicating how rows should be weighted (higher = more rows in resulting matrix).
Raising both alpha and beta equally will find smaller and denser biclusters.
}
  \item{beta}{
The beta term from the BicBin algorithm indicating how columns should be weighted (higher = more columns in resulting matrix).
Raising both alpha and beta equally will find smaller and denser biclusters.
}
  \item{minDensityP1}{
The minimum bicluster density (P_1 in algorithm A_2). Biclustering is run iteratively
until a cluster of this density is found. Setting it to 1 will find non-sparse clusters, and
setting it to 0 will return the best scoring clusters regardless of density.
}
  \item{tries}{
The number of random iterations to try. Users should try multiple values to test for convergence
(i.e. the number needed to produce repeatable results).
}
}
\value{
Returns a list of the bicluster coordinates and scores. The elements
rows and cols contain the row and column numbers for the bicluster, as well as the score for that bicluster.
For example, you can get the coordinates for the first bicluster as follows: bicluster[[1]]$rows, bicluster[[1]]$cols
}
\references{
Biclustering sparse binary genomic data. van Uitert, Meuleman W, Wessels L. J Comput Biol. 2008 Dec;15(10):1329-45. doi: 10.1089/cmb.2008.0066.
}
\author{
Tyler William H Backman
}
\examples{
cores <- 1 # set cpu cores
gc() # clear unused memory
set.seed(123) # get same clusters each time
registerDoMC(cores=cores) # register parallel backend
tries <- cores*10 # make tries an even multiple of cores
row.names(mat) <- as.character(1:nrow(mat))
colnames(mat) <- as.character(1:ncol(mat))
mat <- matrix(sample(c(rep.int(0, 100), 1), 800*250, replace=TRUE), 800, 250)
bicluster <- bicBinCluster(binaryMatrix=mat, totalClusters=1,
	alpha=0.5, beta=0.8, tries=tries, minDensityP1=0)
firstBicluster <- mat[bicluster[[1]]$rows, bicluster[[1]]$cols]
firstBiclusterScore <- bicluster[[1]]$score
}
\keyword{ utilities }
