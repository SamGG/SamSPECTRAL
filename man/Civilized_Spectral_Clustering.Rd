\name{Civilized_Spectral_Clustering}
\Rdversion{1.1}
\alias{Civilized_Spectral_Clustering}

\title{
    Runs the spectral clustering algorithm on the sample points.
}

\description{
    The representatives of communities are considered as the vertices of a graph. Assuming the edges have been weighted according to the equivalent conductance between them, this function runs the classic spectral clustering on the graph.
}

\usage{
Civilized_Spectral_Clustering(full, maximum.number.of.clusters, society, conductance,iterations=200, number.of.clusters=NA,  
								eigenvalues.num =NA, talk=TRUE,stabilizer=1000)
}

\arguments{
  \item{full}{
     The matrix containing the coordinates of all data points. }

  \item{maximum.number.of.clusters}{
    This parameter is used for fitting the regression line.}

  \item{number.of.clusters}{
    The default value is NA which leads to computating the number of spectral clusters automatically, otherwise this number will
		determine the number of spectral clusters. }

  \item{society}{
    The list of communities.}
  
  \item{conductance}{
    A matrix in which each entry is the conductance between two communities.}
  
  \item{iterations}{
  	Number of iterations for the k-means algorithm used by the spectral procedure. 200 is an appropriate value.}

  \item{talk}{
    A boolean flag with default value TRUE. Setting it to FALSE will keep running the procedure quite with no messages.}

  \item{eigenvalues.num}{
        An integer with default value NA which prevents ploting the curve of eigenvalues. Otherwise, they will be ploted upto this number.}

  \item{stabilizer}{
    The larger this integer is, the final results will be more stable because the underlying kmeans will restart many more times.}        
  
}

\value{
    \item{labels.for_num.of.clusters}{
        The k'th element of this list is a vector containing the labels as result of clustering to k parts.}
    \item{number.of.clusters }{
        A list containing the desired cluster numbers.}
    \item{eigen.space}{
        The eigen vectors and eigen values of the normalized adjacency matrix computed for spectral clustering.}
}

\references{
Zare, H. and Shooshtari, P.  and Gupta, A. and Brinkman R.B. (2009). Data Reduction for Spectral Clustering to Analyse High Throughput Flow Cytometry Data. submitted to BMC Bioinformatics. }


\author{
Parisa Shooshtari and Habil Zare}

\seealso{
\code{\link{SamSPECTRAL}}
}

\examples{

	 \dontrun{
	    library(SamSPECTRAL)

	   # Reading data file which has been transformed using log transform
	    data(small_data)
		full <- small
	
		# Parameters:
		m <- 3000; ns <- 200; sl <- 3; cwt <-1; precision <- 6; mnc <-30 
	
    	# Sample the data and build the communities
    	society <- Building_Communities(full=full,m=m, space.length=sl, community.weakness.threshold=cwt)
    	
	
    	# Compute conductance between communities
    	conductance <- Conductance_Calculation(full=full, normal.sigma=ns, space.length=sl, society=society, precision=precision)
    	
    	# Use spectral clustering to cluster the data
		# First example:
	    clust_result <- Civilized_Spectral_Clustering(full=full, maximum.number.of.clusters=mnc, society=society, conductance=conductance)    
    	number.of.clusters <- clust_result@number.of.clusters
    	labels.for_num.of.clusters <- clust_result@labels.for_num.of.clusters
		L <- labels.for_num.of.clusters[[number.of.clusters]]
    	# plot(full, pch='.', col= L)
		

		# Second example:
		number.of.clusters <- c(35,20)	
			# This is faster than runnig Civilized_Spectral_Clustering() twice because the eigen space is not needed to be computed again.     
		clust_result.not.automatic <- 
			Civilized_Spectral_Clustering(full=full, society=society, conductance=conductance, number.of.clusters =number.of.clusters)    
	    labels.for_num.of.clusters <- clust_result.not.automatic@labels.for_num.of.clusters
		L35 <- labels.for_num.of.clusters[[35]]
		L20 <- labels.for_num.of.clusters[[20]]
	    # plot(full, pch='.', col= L35)
	}

}


\keyword{cluster}
\keyword{graphs}

