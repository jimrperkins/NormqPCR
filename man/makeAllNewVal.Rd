\name{makeAllNewVal}
\alias{makeAllNewVal}
\alias{makeAllNewVal,qPCRBatch-method}
\title{ Make all Cq values NA  }
\description{ Make all Cq values for a given detector NA when the number of NAs for that detector is above a given threshold
}
\usage{
makeAllNewVal(qPCRBatch, \dots)

\S4method{makeAllNewVal}{qPCRBatch}(qPCRBatch, contrastM, sampleMaxM, newVal)
}
\arguments{
  \item{qPCRBatch}{ Expression set containing qPCR data.}
  \item{\dots}{ Extra arguments, detailed below}
  \item{contrastM}{ Contrast Matrix like that used in \code{limma}. Columns represent the different samples types, rows are the different samples,
with a 1 or 0 in the matrix indicating which sample types the different samples belong to.}
  \item{sampleMaxM}{ Sample Max Matrix. Columns represent the different sample types. There is one value per column, which represents the max number of NAs allowed for that sample type.}
  \item{newVal}{ New value to give the values in the group where the NAs are above the threshold.}
}
\details{
  Make all a given value when number of NAs above a given threshold, with different thresholds for the different sample classes, using sMaxM and contM to provide this information, as detailed below.
}
\value{
  \code{qPCRBatch} object with a new exprs slot, everything else equal
}
%\references{ }
\author{ James Perkins \email{jperkins@biochem.ucl.ac.uk}}
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
  qPCRBatch.taqman.replaced <- makeAllNewVal(qPCRBatch.taqman, contM, sMaxM)
  exprs(qPCRBatch.taqman.replaced)["Ccl20.Rn00570287_m1",]
}
\keyword{data}
