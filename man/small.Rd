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
  The is an FCS file.}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

\references{
Faithful Sampling for Spectral Clustering to Analyse High Throughput Biological Data.}


\examples{
data(small_data)
	full <- small

    plot(full, pch='.')

}
\keyword{datasets}
