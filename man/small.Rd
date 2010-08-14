\name{small}
\Rdversion{1.1}
\alias{small}
\docType{data}
\title{
Flow cytometry data to test SamSPECTRAL algorithm.}

\description{
This FCS file is a small one used to show how to set SamSPECTRAL parameters.}


\usage{data(small_data)}
\format{
  This is an FCS file.}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

\references{
Zare, H. and Shooshtari, P.  and Gupta, A. and Brinkman R.B. (2009). Data Reduction for Spectral Clustering to Analyse High Throughput Flow Cytometry Data. BMC Bioinformatics, 2010, 11:403. }


\examples{
data(small_data)
	full <- small

    plot(full, pch='.')

}
\keyword{datasets}
