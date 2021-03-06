\name{SamSPECTRAL-package}
\Rdversion{1.1}
\alias{SamSPECTRAL-package}
%\alias{SamSPECTRAL}
\docType{package}

\title{
    Identifying cell populations in flow cytometry data.
}

\description{
    Using a faithful sampling procedure, SamSPECTRAL reduces the size of data points such that applying spectral clustering algorithm on large data such as flow cytometry is possible. Before running the spectral clustering algorithm, it uses potential theory to define similarity between sampled points.
}

\details{
\tabular{ll}{
Package: \tab SamSPECTRAL\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-08-31\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The main function is SamSPECTRAL. It can be loaded using the command library(SamSPECTRAL) in R.
Some parameters should be set properly including: dimensions, normal.sigma and separation.factor. 
These parameters can be adjusted for a data set by running the algorithm on some samples of that data set. (Normally, 2 or 3 samples are sufficient).
Then the function SamSPECTRAL() can be applied to all samples in the data set to identify cell populations in each sample data.
}

\author{
Habil Zare and Parisa Shooshtari 

Maintainer: Habil Zare <hzare@bccrc.ca>
}

\references{
Zare, H. and Shooshtari, P.  and Gupta, A. and Brinkman R.B: Data Reduction for Spectral Clustering to Analyse High Throughput Flow Cytometry Data. BMC Bioinformatics, 2010, 11:403.}



\keyword{cluster}
\seealso{
\code{\link{SamSPECTRAL}}, \code{\link{Building_Communities}}, \code{\link{Conductance_Calculation}}, 
\code{\link{Civilized_Spectral_Clustering}}, \code{\link{Connecting}}
}

\examples{

	\dontrun{
	    library(SamSPECTRAL)

	   # Reading data file which has been transformed using log transform
	    data(small_data)
		full <- small
		
	    L <- SamSPECTRAL(data.points=full,dimensions=c(1,2,3), normal.sigma = 200, separation.factor = 0.39)
	    
	    plot(full, pch='.', col= L)
	}

}

