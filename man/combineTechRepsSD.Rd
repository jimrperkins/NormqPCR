\name{combineTechRepsWithSD}
\alias{combineTechRepsWithSD}
\alias{combineTechRepsWithSD,qPCRBatch-method}
\title{ Combines Technical Replicates }
\description{
Takes expression set of qPCR values containing technical replicates and combines them.
In addition the appropriate standard deviation (SD) is computed.
}
\usage{
combineTechRepsWithSD(qPCRBatch, \dots)

\S4method{combineTechRepsWithSD}{qPCRBatch}(qPCRBatch, calc="arith") 
}
\arguments{
  \item{qPCRBatch}{ Expression set containing qPCR data, read in by a ReadqPCR function and containing technical reps, denoted by \code{_TechRep.n} suffix. }
  \item{\dots}{ Extra arguments, detailed below }
  \item{calc}{ use median, arithmetic or geometric mean for combining the values}
}
\details{
  Takes \code{exprs} of qPCR values containing technical replicates and combines them using a specified centrality measure.
  
  The arithmetic mean (\code{calc="arith"}) is combined with the classical standard deviation.  
  In case of the geometric mean (\code{calc="geom"}) the classical standard deviation of the log-values is exponentiated.
  The median (\code{calc="median"}) is calculated in connection with the MAD.
}
\value{
 \code{qPCRBatch} with same number of samples, but with less features, since all technical replicates are replaced with a single value of their means.
 In addition the slot \code{assayData} includes a matrix with SD values which can be accessed via \code{se.exprs}.
}

\references{ 
  Perkins, JR, Dawes, JM, McMahon, SB, Bennett, DL, Orengo, C, Kohl, M (2012). 
  ReadqPCR and NormqPCR: R packages for the reading, quality checking and 
  normalisation of RT-qPCR quantification cycle (Cq) data. 
  \emph{BMC Genomics}, \bold{13}, 1:296.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{combineTechReps}}}
\examples{
  path <- system.file("exData", package = "NormqPCR")
  qPCR.example.techReps <- file.path(path, "qPCR.techReps.txt")
  qPCRBatch.qPCR.techReps <- read.qPCR(qPCR.example.techReps)
  rownames(exprs(qPCRBatch.qPCR.techReps))
  combinedTechReps <- combineTechRepsWithSD(qPCRBatch.qPCR.techReps)
  rownames(exprs(combinedTechReps))
  exprs(combinedTechReps)
  se.exprs(combinedTechReps)
}
\keyword{data}
