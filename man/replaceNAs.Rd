\name{replaceNAs}
\alias{replaceNAs}
\alias{replaceNAs,qPCRBatch-method}
\title{ Replace NAs with a given value }
\description{ Replace NAs with a given value
}
\usage{
replaceNAs(qPCRBatch, \dots)

\S4method{replaceNAs}{qPCRBatch}(qPCRBatch, newNA)
}
\arguments{
  \item{qPCRBatch}{ Expression set containing qPCR data.}
  \item{\dots}{ Extra arguments, detailed below}
  \item{newNA}{ The new value to replace the NAs with}
}
\details{
  Replaces NA values in the exprs slot of the qPCRBatch object with a given number
}
\value{
  qPCRBatch object with a new exprs slot
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
  qPCRBatch.taqman.replaced <- replaceNAs(qPCRBatch.taqman, newNA = 40)
  exprs(qPCRBatch.taqman.replaced)["Ccl20.Rn00570287_m1",]
}
\keyword{data}
