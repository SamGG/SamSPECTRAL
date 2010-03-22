\name{Connecting}
\Rdversion{1.1}
\alias{Connecting}

\title{
    Combines the spectral clusters to build the connected components.
}

\description{
    Considering some biological criterion based on density, the clusters which are identified by spectral clustering are combined to estimate biological populations.}

\usage{
Connecting(full, society,conductance, number.of.clusters, labels.for_num.of.clusters, separation.factor, talk=TRUE)
}


\arguments{
  \item{full}{
    The matrix containing the coordinates of all data points. }
  
  \item{society}{  
    The list of communities.}
    
  \item{conductance}{
    A matrix in which each entry is the conductance between two communities.}

  \item{number.of.clusters}{
    A list containing the desired cluster numbers.}
  
  \item{labels.for_num.of.clusters}{
    The k'th element of this list, is a vector containing the labels as result of clustering to k parts.}

  \item{separation.factor}{
    This threshold controls to what extend clusters should be combined or kept separate.}   

  \item{talk}{
    A boolean flag with default value TRUE. Setting it to FALSE will keep running the procedure quite with no messages.}
}

\details{
	A hint for setting \code{separation.factor}: While \code{separation.factor}=0.7 is normally an appropriate value for many datasets, 
		for others some value in range 0.3 to 1.2 may produce better results depending on what populations are of particular interest.
}

\value{
    Returns a vector containing the labels that determines to which component each data point belongs.
}

\references{
Faithful Sampling for Spectral Clustering to Analyse High Throughput Biological Data.}

\author{
Parisa Shooshtari and Habil Zare 
}

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
	    clust_result <- Civilized_Spectral_Clustering(full=full, maximum.number.of.clusters=mnc, society=society, conductance=conductance)    
    	number.of.clusters <- clust_result@number.of.clusters
    	labels.for_num.of.clusters <- clust_result@labels.for_num.of.clusters
		L <- labels.for_num.of.clusters[[number.of.clusters]]
    	# plot(full, pch='.', col= L)
		
    	
    	# Connect components
    	L <- Connecting(full=full, society=society, conductance=conductance, number.of.clusters=number.of.clusters,
			 labels.for_num.of.clusters=labels.for_num.of.clusters, separation.factor=0.39)
    	
    	plot(full, pch='.', col= L)

	}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}


