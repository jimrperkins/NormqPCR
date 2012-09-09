\name{BladderRepro}
\alias{BladderRepro}
\alias{BladderRepro.qPCRBatch}
\docType{data}
\title{ Dataset of Andersen et al (2004) }
\description{
  This dataset was used in Andersen et al (2004) to demonstrate normalization 
  of real-time quantitative RT-PCR data by geometric averaging of housekeeping genes.
}
\usage{data(BladderRepro)}
\format{
  A qPCRBatch object which contains an expression matrix with the expression of
  8 genes measured in 26 samples. The sample information is saved in the
  \code{phenoData} slot with variables
  \describe{
    \item{\code{Sample.no.}}{sample number. }
    \item{\code{Grade}}{Grade of bladder cancer. }
  }
  The following information on the measured genes is saved in the variables 
  \code{Symbol} and \code{Gene.name} of the \code{featureData} slot.
  \describe{
    \item{\code{CD14}}{CD14 antigen. }
    \item{\code{FCN1}}{Ficolin (collagen/fibrinogen domain containing) 1. }
    \item{\code{CCNG2}}{Cyclin G2. }
    \item{\code{NPAS2}}{Neuronal PAS domain protein 2. }
    \item{\code{UBC}}{Ubiquitin C. }
    \item{\code{CFL1}}{Cofilin 1 (non-muscle). }
    \item{\code{ACTB}}{Actin, beta. }
    \item{\code{GAPD}}{Glyceraldehyde-3-phosphate dehydrogenase. }
  }
  For a detailed annotation see Table 1 and Supplementary table 1 in Anderson et al. (2004).
}
\details{
  This data set was used to check the reproducibility of the results obtained
  in Andersen et al (2004).
}
\source{
  The data set was obtained from
  \url{http://www.mdl.dk/Publications_sup1.htm}
}
\references{
  Claus Lindbjerg Andersen, Jens Ledet Jensen and Torben Falck Orntoft (2004).
  Normalization of Real-Time Quantitative Reverse Transcription-PCR Data: 
  A Model-Based Variance Estimation Approach to Identify Genes Suited for
  Normalization, Applied to Bladder and Colon Cancer Data Sets. 
  CANCER RESEARCH 64, 5245-5250, August 1, 2004.
  \url{http://cancerres.aacrjournals.org/cgi/content/full/64/15/5245}
}
\examples{
  data(BladderRepro)
  BladderRepro
  head(exprs(BladderRepro))
  pData(BladderRepro)
  fData(BladderRepro)
}
\keyword{datasets}
