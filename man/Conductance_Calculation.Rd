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
Conductance_Calculation(full, normal.sigma, space.length, society, precision, talk=TRUE)
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
}

\value{
    Returns a matrix in which each entry is the conductance between two communities.}
    
\references{
Faithful Sampling for Spectral Clustering to Analyse High Throughput Biological Data.}

\author{
Parisa Shooshtari and Habil Zare}


\seealso{
\code{\link{SamSPECTRAL}}
}
\examples{

	\dontrun{
	    library(SamSPECTRAL)
    	data(stm)
	   # Read data files and transform them using log transform
		data.points <- exprs(stmFSC) 
		dimensions <- c(3,4,7)
		full <- log10(data.points[,dimensions])
	
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

