\name{deltaCt}
\alias{deltaCt}
\alias{deltaCq}
\alias{deltaCt-methods}
\alias{deltaCt,qPCRBatch-method}
\title{ Perform normalization with a given housekeeping gene }
\description{
Normalise qPCR eset using a given housekeeping gene as control, then perform differential expression analysis using the delta delta Ct method
}
\usage{
deltaCt(qPCRBatch, \dots)
\S4method{deltaCt}{qPCRBatch}(qPCRBatch, hkgs, combineHkgs=FALSE, calc="arith")
deltaCq(qPCRBatch, hkgs, combineHkgs, calc)
}
\arguments{
  \item{qPCRBatch}{ qPCR-specific expression set, containing qPCR data. }
  \item{\dots}{ Extra arguments, detailed below }
  \item{hkgs}{ String containing the name of the name of the housekeeping gene which will be used to normalise the rest of the genes.}
  \item{combineHkgs}{  Logical - if TRUE, then as long as more than one housekeeper given for argument hkgs, it will combine the housekeepers by finding the geometric mean. Housekeepers can be found using geNorm or NormFinder algorithms.}
  \item{calc}{ use arithmetic or geometric mean.}
}
\details{
  Takes expression set of qPCR values and normalises them using a housekeeping gene. Returns a qPCRBatch with exprs set of the same dimensions but with the given hkg value subtracted.
}
\value{
 qPCRBatch with exprs set of the same dimensions but with the given hkg value subtracted.

}
\references{
  Kenneth Livak, Thomase Schmittgen (2001).
  Analysis of Relative Gene Expression Data Using Real-Time Quantitative
  PCR and the 2^DDCt Method.
  Methods 25, 402-408, 2001
  \url{http://www.ncbi.nlm.nih.gov/pubmed/11846609}
}
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk}}
%\note{}
\seealso{\code{selectHKs}, \code{deltaDeltaCq}}
\examples{
  path <- system.file("exData", package = "NormqPCR")
  taqman.example <- file.path(path, "example.txt")
  qPCRBatch.taqman <- read.taqman(taqman.example)
  hkgs<-"Actb-Rn00667869_m1"
  qPCRBatch.norm <- deltaCq(qPCRBatch =  qPCRBatch.taqman, hkgs = hkgs, calc="arith")
  head(exprs(qPCRBatch.norm))
}
\keyword{data}
