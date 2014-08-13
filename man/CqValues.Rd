\name{CqValues}
\alias{CqValues}
\alias{CqValues-methods}
\alias{CqValues,CyclesSet-method}
\title{ Compute Cq value and amplification efficiency }
\description{
  This function calculates Cq value and amplification efficiency for a \code{CyclesSet}. 
  It is based on function \code{\link[qpcR]{pcrbatch}} of package \pkg{qpcR}. 
}
\usage{
CqValues(object, ...)
\S4method{CqValues}{CyclesSet}(object, Effmethod = "expfit", group = NULL, 
                      model = l5, check = "uni2", checkPAR = parKOD(), 
                      remove = "none", exclude = NULL, type = "cpD2", 
                      labels = NULL, norm = FALSE, baseline = FALSE, 
                      basefac = 1, smooth = NULL, 
                      smoothPAR = list(span = 0.1), 
                      factor = 1, opt = FALSE, 
                      optPAR = list(sig.level = 0.05, crit = "ftest"), 
                      plot = FALSE, verbose = FALSE, \dots)
}
\arguments{
  \item{object}{an object of class \code{\link[ReadqPCR]{CyclesSet}}.}
  \item{Effmethod}{a character vector defining the methods for computing amplification efficiency.}
  \item{group}{a vector containing the grouping for possible replicates.}
  \item{model}{the model to be used for all runs. Default model is \code{l5}.}
  \item{check}{the method for kinetic outlier detection in \code{\link[qpcR]{KOD}}. Method \code{"uni2"} is set as default which is a test on sigmoidal structure.}
  \item{checkPAR}{parameters to be supplied to the \code{check} method. See \code{\link[qpcR]{parKOD}}.}
  \item{remove}{indicates which runs to be removed. Either \code{none} of them, those which failed to \code{fit} or from the outlier methods.}
  \item{exclude}{indicates samples to be excluded from calculation, either "" for samples with missing column names or a regular expression defining columns (samples); 
                 see 'Details' and 'Examples' in \code{\link[qpcR]{modlist}}.}
  \item{type}{the point on the amplification curve which is used for efficiency estimation; see \code{\link[qpcR]{efficiency}}.}
  \item{labels}{a vector containing labels which define replicate groups. See more details in \code{\link[qpcR]{pcrbatch}} and \code{\link[qpcR]{ratiobatch}}.}
  \item{norm}{a logical value which determines whether the raw data should be normalized within [0, 1] before model fitting or not.}
  \item{baseline}{type of baseline subtraction. More details in \code{\link[qpcR]{efficiency}}.}
  \item{basefac}{a factor when using averaged baseline cycles, such as \code{0.95}.}
  \item{smooth}{the curve smoothing method. See more details in \code{\link[qpcR]{pcrbatch}}.}
  \item{smoothPAR}{parameters to be supplied to smoothing method in \code{smooth}.}
  \item{factor}{a multiplication factor for the fluorescence response values.}
  \item{opt}{a logical value which determines whether model selection should be applied to each model or not.}
  \item{optPAR}{parameters to be supplied for model selection in \code{\link[qpcR]{mselect}}.}
  \item{plot}{a logical value. If \code{TRUE}, the single runs are plotted from the internal \code{modlist} for diagnostics.}
  \item{verbose}{a logical value. If \code{TRUE}, fitting and tagging results will be displayed in the console.}
  \item{\dots}{other parameters to be passed to downstream methods.}
}
\details{
  Allows the user to compute Cq value and amplification efficiency. In addition, all values generated during the computations are saved. 
  This function has four choices of methods for computing amplification efficiency values which are the methods provided by package \pkg{qpcR}. 
  
  More details on technical replication and normalization is given in the vignette \code{NormqPCR}.
}
\value{Object of class \code{"qPCRBatch"}.}
\references{
  Perkins, JR, Dawes, JM, McMahon, SB, Bennett, DL, Orengo, C, Kohl, M (2012). 
  ReadqPCR and NormqPCR: R packages for the reading, quality checking and 
  normalisation of RT-qPCR quantification cycle (Cq) data. 
  \emph{BMC Genomics}, \bold{13}, 1:296.
}
\author{ Nor Izayu Abdul Rahman, Matthias Kohl \email{Matthias.Kohl@stamats.de} }
%\note{
%}
\seealso{ \code{\link[qpcR]{pcrbatch}}, \code{\link[ReadqPCR]{CyclesSet-class}}, \code{\link[ReadqPCR]{qPCRBatch-class}}}
\examples{
  ## Read in the raw qPCR data from file "LC480_Example.txt"
  path <- system.file("exData", package = "ReadqPCR")
  LC480.example <- file.path(path, "LC480_Example.txt")
  cycData <- read.LC480(file = LC480.example)

  ## Read in the sample information data from file "LC480_Example_SampleInfo.txt". 
  LC480.SamInfo <- file.path(path, "LC480_Example_SampleInfo.txt")
  samInfo <- read.LC480SampleInfo(LC480.SamInfo)

  ## Merge information
  cycData1 <- merge(cycData, samInfo)

  ## Compute Cq values
  ## 1) use sigmoidal model
  res1 <- CqValues(cycData1, Effmethod = "sigfit")
  res1
  effs(res1)
  se.effs(res1)

  ## 2) fit exponential model (default)
  res2 <- CqValues(cycData1, Effmethod = "expfit")
  res2
  effs(res2)
  se.effs(res2)

  ## 3) use window of linearity
  res3 <- CqValues(cycData1, Effmethod = "sliwin")
  res3
  effs(res3)
  se.effs(res3)

  ## 4) linear regression of efficiency
  res4 <- CqValues(cycData1, Effmethod = "LRE")
  res4
  effs(res4)
  se.effs(res4)
}
\keyword{classes}
