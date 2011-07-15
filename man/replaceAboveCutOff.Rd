\name{replaceAboveCutOff}
\alias{replaceAboveCutOff}
\alias{replaceAboveCutOff,qPCRBatch-method}
\title{ Replace Ct values with new value}
\description{ Replace Ct values above a given threshold with a new value
}
\usage{
replaceAboveCutOff(qPCRBatch, \dots)

\S4method{replaceAboveCutOff}{qPCRBatch}(qPCRBatch, newVal=NA, cutOff=38)
}
\arguments{
  \item{qPCRBatch}{ Expression set containing qPCR data. }
  \item{\dots}{ Extra arguments, detailed below }
  \item{newVal}{ The new value with which to replace the values above the cutoff }
  \item{cutOff}{ the minimal threshold above which the values will be replaced }
}
\details{
  Replaces values in the exprs slot of the \code{qPCRBatch} object that are above a threshold value with a new number
}
\value{
  \code{qPCRBatch} object with a new exprs slot
}
%\references{}
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk}}
%\note{}
%\seealso{}
\examples{
  path <- system.file("exData", package = "NormqPCR")
  taqman.example <- file.path(path, "example.txt")
  qPCRBatch.taqman <- read.taqman(taqman.example)
  exprs(qPCRBatch.taqman)["Ccl20.Rn00570287_m1",]
  qPCRBatch.taqman.replaced <- replaceAboveCutOff(qPCRBatch.taqman, newVal = NA, cutOff = 35)
  exprs(qPCRBatch.taqman.replaced)["Ccl20.Rn00570287_m1",]
}
\keyword{data}
