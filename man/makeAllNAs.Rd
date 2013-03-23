\name{makeAllNAs}
\alias{makeAllNAs}
\alias{makeAllNAs,qPCRBatch-method}
\title{ Make all Cq values NA  }
\description{ Make all Cq values for a given detector NA when the number of NAs for that detector is above a given threshold
}
\usage{
makeAllNAs(qPCRBatch, \dots)

\S4method{makeAllNAs}{qPCRBatch}(qPCRBatch, contrastM, sampleMaxM)
}
\arguments{
  \item{qPCRBatch}{ Expression set containing qPCR data.}
  \item{\dots}{ Extra arguments, detailed below}
  \item{contrastM}{ Contrast Matrix like that used in \code{limma}. Columns represent the different samples types, rows are the different samples,
with a 1 or 0 in the matrix indicating which sample types the different samples belong to.}
  \item{sampleMaxM}{ Sample Max Matrix. Columns represent the different sample types. There is one value per column, which represents the max number of NAs allowed for that sample type.}
}
\details{
  Make all NAs when number of NAs above a given threshold
}
\value{
  \code{qPCRBatch} object with a new exprs slot, everything else equal
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
  # read in the data
  path <- system.file("exData", package = "NormqPCR")
  taqman.example <- file.path(path, "example.txt")
  qPCRBatch.taqman <- read.taqman(taqman.example)
  exprs(qPCRBatch.taqman)["Ccl20.Rn00570287_m1",] # values before

  # make contrastM
  a <- c(0,0,1,1,0,0,1,1) # one for each sample type, with 1 representing
  b <- c(1,1,0,0,1,1,0,0) # position of sample type in the samplenames vector
  contM <- cbind(a,b)
  colnames(contM) <- c("case","control") # then give the names of each sample type
  rownames(contM) <- sampleNames(qPCRBatch.taqman) # and the rows of the matrix
  contM

  # make sampleMaxM
  sMaxM <- t(as.matrix(c(3,3))) # now make the sample max matrix
  colnames(sMaxM) <- c("case","control") # make sure these line up with samples
  sMaxM

  # function
  qPCRBatch.taqman.replaced <- makeAllNAs(qPCRBatch.taqman, contM, sMaxM)
  exprs(qPCRBatch.taqman.replaced)["Ccl20.Rn00570287_m1",]
}
\keyword{data}
