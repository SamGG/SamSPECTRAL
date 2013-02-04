\name{kneepointDetection}
\Rdversion{1.1}
\alias{kneepointDetection}

\title{
  Fits 2 regression lines to data to estimate the knee (or elbow) point.
}

\description{
  With an appropriate sigma value, the curve of eigenvalues has a knee point shape. The bending point is a good estimate for the number of informative spectral clusters because the eigenvalues above the corresponding threshold can reasonably be assumed to be close to 1. This function estimate the knee point by fitting 2 lines using linear regression.}

\usage{
  kneepointDetection(vect, PlotFlag=FALSE)
}


\arguments{
  \item{vect}{
    The vector of values on which the 2 regression lines will be fitted.}
  
  \item{PlotFlag}{  
    If TRUE and in unix, an animation will be produced in tmpfigs folder 
      that shows how the best selected model in gif format.}
}

\details{
  The running time is in order of minutes for 100 points. This function was borrowed from flowMeans package and for application in SamSPECTRAL package, it was customized such that the first line is always horizontal. 
}

\value{
  Returns a list where MinIndex is the index of the knee point and l1 and l2 the fitted lines.
}

\references{
  Aghaeepour N., Nikolic R., Hoos HH., Brinkman RR.: Rapid cell population identification in flow cytometry data. Cytometry A, 2011, 79:6.}



\author{Nima Aghaeepour}

\seealso{
\code{\link[flowMeans]{changepointDetection}}
}

\examples{
  ## Data
    values <- rep(1,times=10)
    values <- c(values,(10:0)/10)
    
  ## Looks like knee point:
  plot(values)

  ## Find the knee point:
  detected <- kneepointDetection(vect=values, PlotFlag=FALSE)
  print(detected$MinIndex)
  ##  Also, under unix, set PlotFlag=TRUE and look at animation.gif.
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}


