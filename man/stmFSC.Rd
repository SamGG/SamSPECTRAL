\name{stmFSC}
\Rdversion{1.1}
\alias{stmFSC}
\docType{data}
\title{
Flow cytometry data to test SamSPECTRAL algorithm.}

\description{
This FCS file is used as demo data to illustrate SamSPECTRAL capabilities in identifying cell populations.}


\usage{data(stm)}
\format{
  The is an FCS file.}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

\references{
Zare, H. and Shooshtari, P.  and Gupta, A. and Brinkman R.B. (2009). Data Reduction for Spectral Clustering to Analyse High Throughput Flow Cytometry Data. BMC Bioinformatics, 2010, 11:403. }

\examples{
data(stm)
   # Read data files and transform them using log transform
	data.points <- stmFSC@exprs  
	dimensions <- c(3,4,7)
	full <- log10(data.points[,dimensions])

    plot(full, pch='.')

}
\keyword{datasets}
