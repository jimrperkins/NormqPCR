\name{replaceAboveCutOff}
\alias{replaceAboveCutOff}
\alias{replaceAboveCutOff,qPCRBatch-method}
\title{ Replace Cq values with new value}
\description{ Replace Cq values above a given threshold with a new value
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
\references{ 
  Perkins, JR, Dawes, JM, McMahon, SB, Bennett, DL, Orengo, C, Kohl, M (2012). 
  ReadqPCR and NormqPCR: R packages for the reading, quality checking and 
  normalisation of RT-qPCR quantification cycle (Cq) data. 
  \emph{BMC Genomics}, \bold{13}, 1:296.
}
\author{ James Perkins \email{jimrperkins@gmail.com}}
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
