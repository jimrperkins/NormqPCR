\name{combineTechReps}
\alias{combineTechReps}
\alias{combineTechReps,qPCRBatch-method}
\title{ Combines Technical Replicates }
\description{
Takes expression set of qPCR values containing technical replicates and combines them.
}
\usage{
combineTechReps(qPCRBatch, \dots)

\S4method{combineTechReps}{qPCRBatch}(qPCRBatch, calc="arith") 
}
\arguments{
  \item{qPCRBatch}{ Expression set containing qPCR data, read in by a ReadqPCR function and containing technical reps, denoted by \code{_TechRep.n} suffix. }
  \item{\dots}{ Extra arguments, detailed below }
  \item{calc}{ use median, arithmetic or geometric mean for combining the values}
}
\details{
  Takes \code{exprs} of qPCR values containing technical replicates and combines them using a specified centrality measure.
}
\value{
 \code{qPCRBatch} with same number of samples, but with less features, since all technical replicates are replaced with a single value of their means.
}

\references{ 
  Perkins, JR, Dawes, JM, McMahon, SB, Bennett, DL, Orengo, C, Kohl, M (2012). 
  ReadqPCR and NormqPCR: R packages for the reading, quality checking and 
  normalisation of RT-qPCR quantification cycle (Cq) data. 
  \emph{BMC Genomics}, \bold{13}, 1:296.
}
\author{ James Perkins \email{j.perkins@ucl.ac.uk}}
%\note{}
%\seealso{}
\examples{
  path <- system.file("exData", package = "NormqPCR")
  qPCR.example.techReps <- file.path(path, "qPCR.techReps.txt")
  qPCRBatch.qPCR.techReps <- read.qPCR(qPCR.example.techReps)
  rownames(exprs(qPCRBatch.qPCR.techReps))
  combinedTechReps <- combineTechReps(qPCRBatch.qPCR.techReps)
  rownames(exprs(combinedTechReps))
}
\keyword{data}
