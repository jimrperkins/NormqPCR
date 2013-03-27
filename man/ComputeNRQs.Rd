\name{ComputeNRQs}
\alias{ComputeNRQs}
\alias{ComputeNRQs-methods}
\alias{ComputeNRQs,qPCRBatch-method}
\title{ Compute Normalized Relative Quantities }
\description{
  This function computes normalized relative quantities (NRQs) for a \code{qPCRBatch}. 
}
\usage{
ComputeNRQs(qPCRBatch, ...)
\S4method{ComputeNRQs}{qPCRBatch}(qPCRBatch, hkgs)
}
\arguments{
  \item{qPCRBatch}{an object of class \code{\link[ReadqPCR]{qPCRBatch}}.}
  \item{hkgs}{Names of reference/housekeeping genes.}
  \item{\dots}{other parameters to be passed to downstream methods.}
}
\details{
  Allows the user to normalized relative quantities as defined in Hellemanns et al. (2007).
}
\value{Object of class \code{"qPCRBatch"}.}
\references{
  Jan Hellemans, Geert Mortier, Anne De Paepe, Frank Speleman and Jo Vandesompele (2007).
  qBase relative quantification framework and software for management and automated analysis 
  of real-time quantitative PCR data.
  \emph{Genome Biology}, 8:R19

  Perkins, JR, Dawes, JM, McMahon, SB, Bennett, DL, Orengo, C, Kohl, M (2012). 
  ReadqPCR and NormqPCR: R packages for the reading, quality checking and 
  normalisation of RT-qPCR quantification cycle (Cq) data. 
  \emph{BMC Genomics}, \bold{13}, 1:296.
}
\author{ Nor Izayu Abdul Rahman, Matthias Kohl \email{Matthias.Kohl@stamats.de} }
%\note{
%}
\seealso{ \code{\link[ReadqPCR]{qPCRBatch-class}}}
\examples{
  ## Example data
  path <- system.file("exData", package = "ReadqPCR")
  qPCR.example <- file.path(path, "qPCR.example.txt")
  Cq.data <- read.qPCR(qPCR.example)

  ## combine technichal replicates
  Cq.data1 <- combineTechRepsWithSD(Cq.data)

  ## add efficiencies
  Effs <- file.path(path, "Efficiencies.txt")
  Cq.effs <- read.table(file = Effs, row.names = 1, header = TRUE)
  rownames(Cq.effs) <- featureNames(Cq.data1)
  effs(Cq.data1) <- as.matrix(Cq.effs[,"efficiency",drop = FALSE])
  se.effs(Cq.data1) <- as.matrix(Cq.effs[,"SD.efficiency",drop = FALSE])
  
  ##
  res <- ComputeNRQs(Cq.data1, hkgs = c("gene_az", "gene_gx"))
  ## NRQs
  exprs(res)
  ## SD of NRQs
  se.exprs(res)
}
\keyword{classes}
