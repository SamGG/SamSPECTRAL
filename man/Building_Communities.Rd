\name{Building_Communities}
\Rdversion{1.1}
\alias{Building_Communities}

\title{
    Builds the communities from the set of all data points.
}

\description{
    Some sample points are picked up and the points close to each sample point are considered as members of that community.    
}

\usage{
Building_Communities(full, m=3000, space.length=1, community.weakness.threshold=1, talk=TRUE, do.sampling=TRUE,
			   replace.inf.with.extremum=TRUE)
}

\arguments{
  \item{full}{
    The matrix containing the coordinates of all data points. }
    
  \item{m}{
    An integer determining upper and lower bounds on the final number of sample points which will be in range  \code{.95*m/2} and 2 \code{1.1*m}}
    
  \item{space.length}{
    An estimate for the length of a cube that is assumed to contain all data points.}
    
  \item{community.weakness.threshold}{
    The communities with number of members less than this threshold will be ignored. Normally, setting it to 1 is reasonable.}

  \item{talk}{
    A boolean flag with default value TRUE. Setting it to FALSE will keep running the procedure quite with no messages.}

  \item{do.sampling}{
    A boolean flag with default value TRUE. If set to FALSE, the sampling stage will be ignored by picking up all the data points.}

  \item{replace.inf.with.extremum}{
    If TRUE, the Inf and -Inf values will be replaced by maximum and minimum of data in each direction.}
}

\value{
    Returns a society which is a list of communities.}

\references{
Zare, H. and Shooshtari, P.  and Gupta, A. and Brinkman R.B: Data Reduction for Spectral Clustering to Analyse High Throughput Flow Cytometry Data. BMC Bioinformatics, 2010, 11:403.}

\author{
Habil Zare and Parisa Shooshtari 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SamSPECTRAL}},\code{\link{check.SamSPECTRAL.input}}
}

\examples{
	
	 \dontrun{
	    library(SamSPECTRAL)

	   # Reading data file which has been transformed using log transform
	    data(small_data)
		full <- small
	
		# Parameters:
		m <- 3000; ns <- 200; sl <- 3; cwt <-1
	
    	# Sample the data and build the communities
    	society <- Building_Communities(full=full,m=m, space.length=sl, community.weakness.threshold=cwt)
	
		# Ploting the representatives:
		plot(full[society$representatives,])

	}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}


