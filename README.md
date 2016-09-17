# BicBin
Biclustering Sparse Binary Data by Tyler Backman

Installation:
> library(devtools)
> install_github("TylerBackman/BicBin")

This R package is based on the BicBin algorithm published by 
Miranda van Uitert, Wouter Meuleman, and Lodewyk Wessels.
It makes use of the R packages foreach, doMC, and doRNG
to find biclusters in parallel across a large number of CPU cores.

If you use this software please cite the following paper:
Biclustering sparse binary genomic data.
van Uitert, Meuleman W, Wessels L.
J Comput Biol. 2008 Dec;15(10):1329-45. doi: 10.1089/cmb.2008.0066.
http://www.ncbi.nlm.nih.gov/pubmed/19040367

As the original code was published without a license, I am assuming it is in the public
domain. If this is not the case, please contact me immediately and I will take this
software offline.
