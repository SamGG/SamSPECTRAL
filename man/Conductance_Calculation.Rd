\name{Conductance_Calculation}
\Rdversion{1.1}
\alias{Conductance_Calculation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the conductance between communities.
}
\description{
For each two communities, the conductance between their members is summed up and the result is returned as the conductance between the two communities.
}

\usage{
Conductance_Calculation(full, normal.sigma, space.length, society, precision, talk=TRUE, beta=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{full}{
 The matrix containing the coordinates of all data points. }

  \item{normal.sigma}{
    The scaling parameter, the larger it is the algorithm will find smaller clusters.}
    
  \item{space.length}{
    An estimate for the length of a cube that is assumed to contain all data points.}

  \item{society}{
    The list of communities.}
    
  \item{precision}{
    Determines the precision of computations. Setting it to 6 will work and increasing it does not improve results.}

  \item{talk}{
    A boolean flag with default value TRUE. Setting it to FALSE will keep running the procedure quite with no messages.}

  \item{beta}{
    A parameter with default value 4 which must NOT be changed except for huge samples with more than 100,000 data points or for developmental purposes.
 	Setting beta to zero will reduce computational time by applying the following approximation to the conductance calculation step. 
	For each two community, the conductance will be the conductance between their representatives times their sizes.}
}

\value{
    Returns a matrix in which each entry is the conductance between two communities.}
    
\references{
Zare, H. and Shooshtari, P.  and Gupta, A. and Brinkman R.B: Data Reduction for Spectral Clustering to Analyse High Throughput Flow Cytometry Data. BMC Bioinformatics, 2010, 11:403.}



\author{
Habil Zare and Parisa Shooshtari}


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
		m <- 3000; ns <- 200; sl <- 3; cwt <-1; precision <- 6
	
    	# Sample the data and build the communities
    	society <- Building_Communities(full=full,m=m, space.length=sl, community.weakness.threshold=cwt)
    	
	
    	# Compute conductance between communities
    	conductance <- Conductance_Calculation(full=full, normal.sigma=ns, space.length=sl, society=society, precision=precision)
 
	}
    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

