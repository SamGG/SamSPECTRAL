\name{check.SamSPECTRAL.input}
\Rdversion{1.1}
\alias{check.SamSPECTRAL.input}

\title{
    Checks the input to SamSPECTRAL.
}

\description{
    The input to SamSPECTRAL should be a numeric matrix WITHOUT any NA, NaN, and +/- Inf. This function checks the input matrix 
    and prodeuces an error for an inappropriate input.
}

\usage{
check.SamSPECTRAL.input(data.points,dimensions=1:ncol(data.points),replace.inf.with.extremum=FALSE)
}

\arguments{
  \item{data.points}{
	A matrix that contains coordinates of the data points.}

  \item{dimensions}{
    A vector that determines which dimension of the data point matrix are chosen for investigation.}

  \item{replace.inf.with.extremum}{
    If TRUE, the Inf and -Inf values will be replaced by maximum and minimum of data in each direction.}
}


\value{
    Returns a list with the following entries:
    \item{data.matrix}{The data with infinite elements fixed if \code{replace.inf.with.extremum=TRUE}}
    \item{dimensions}{All the checked dimensions.}
    \item{infinite}{Will be TRUE if data contained infinite entries.}
}

\references{
Zare, H. and Shooshtari, P.  and Gupta, A. and Brinkman R.B: Data Reduction for Spectral Clustering to Analyse High Throughput Flow Cytometry Data. BMC Bioinformatics, 2010, 11:403.}



\author{Habil Zare}

\seealso{
\code{\link{SamSPECTRAL}}, \code{\link{Building_Communities}}, \code{\link{Conductance_Calculation}}, 
\code{\link{Connecting}}
}

\examples{

	\dontrun{
 	   library(SamSPECTRAL)
	
 	  # Reading data file which has been transformed using log transform
 	   data(small_data)
		full <- small
		
 	   checked <- check.SamSPECTRAL.input(data.points=full,dimensions=c(1,2,3),replace.inf.with.extremum=TRUE)
 	   
 	   plot(checked$data.matrix, pch='.')
	}    
}

\keyword{cluster}


