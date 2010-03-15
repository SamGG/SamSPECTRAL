\name{SamSPECTRAL}
\Rdversion{1.1}
\alias{SamSPECTRAL}

\title{
    Identifies the cell populations in flow cytometry data.
}

\description{
    Given an FCS file as input, SamSPECTRAL first builds the communities to sample the data points. Then, it builds a graph and after weighting the edges of the graph by conductance computation, it is passed to a classic spectral clustering algorithm to find the spectral clusters. The last stage of SamSPECTRAL is to combine the spectral clusters. The resulting "connected components" estimate biological cell populations in the data sample.
}

\usage{
SamSPECTRAL(data.points, dimensions=colnames(data.points), normal.sigma = 500, separation.factor = 0.8,number.of.clusters = NA, 
	talk = TRUE, precision = 6, eigenvalues.num =NA, return_only.labels=TRUE, do.sampling=TRUE)
}

\arguments{
  \item{data.points}{
	A matrix that contains coordinates of the data points.}

  \item{dimensions}{
    A vector that determines which dimension of the data point matrix are chosen for investigation.}

  \item{normal.sigma}{
	A scaling parameter that determines the "resolution" in the spectral clustering stage. By increasing it, more spectral clusters are identified.
	This can be usefull when "small" population are aimed.} 

  \item{separation.factor}{
    This threshold controls to what extend clusters should be combined or kept separate.}

  \item{number.of.clusters}{
    The default value is NA which leads to computating the number of spectral clusters automatically, otherwise this number will
		determine the number of spectral clusters. }

  \item{talk}{
    A boolean flag with default value TRUE. Setting it to FALSE will keep running the procedure quite with no messages.}

  \item{precision}{
    Determines the precision of computations. Setting it to 6 will work and increasing it does not improve results.}

  \item{eigenvalues.num}{
        An integer with default value NA which prevents ploting the curve of eigenvalues. Otherwise, they will be ploted upto this number.}

  \item{return_only.labels}{
        A boolean flag with default value TRUE. If the user set it to FALSE, SamSPECTRAL function will return all the intermediate objects 
		that are computed during the sampling, similarity calculation, spectral clustering and combining stages.}

  \item{do.sampling}{
    A boolean flag with default value TRUE. If set to FALSE, the sampling stage will be ignored by picking up all the data points.}

}

\details{
	Hints for setting \code{separation.factor} and \code{normal.sigma}: While \code{separation.factor}=0.7 is normally an appropriate value for many datasets, 
		for others some value in range 0.3 to 1.2 may produce better results depending on what populations are of particular interest.
		The larger \code{normal.sigma} is the algorithm will find smaller clusters. It can be adjusted best by considering the plot of eigenvalues as explained in the vignette.
}


\value{
    Returns a vector of labels for data points. If the input parameter return_only.labels is set to FALSE, all the objects 
		that are computed during the intermediate will be returned including: society for sampling stage, conductance for similarity calculation, and
		clustering_result.
}

\references{
Faithful Sampling for Spectral Clustering to Analyse High Throughput Biological Data.}

\author{
Habil Zare and Parisa Shooshtari}

\seealso{
\code{\link{SamSPECTRAL}}, \code{\link{Building_Communities}}, \code{\link{Conductance_Calculation}}, 
\code{\link{Civilized_Spectral_Clustering}}, \code{\link{Connecting}}
}

\examples{

    library(SamSPECTRAL)

   # Reading data file which has been transformed using log transform
    data(small_data)
	full <- small
	
    L <- SamSPECTRAL(data.points=full,dimensions=c(1,2,3), normal.sigma = 200, separation.factor = 0.39)
    
    plot(full, pch='.', col= L)
}

\keyword{cluster}


