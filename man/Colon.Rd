\name{Colon}
\alias{Colon}
\docType{data}
\title{ Colon dataset of Andersen et al (2004) }
\description{
  This dataset was used in Andersen et al (2004) to demonstrate normalization 
  of real-time quantitative RT-PCR data by geometric averaging of housekeeping genes.
}
\usage{data(Colon)}
\format{
  A qPCRBatch object which contains an expression matrix with the expression of
  13 genes measured in 40 samples. The sample information is saved in the
  \code{phenoData} slot with variables
  \describe{
    \item{\code{Sample.no.}}{sample number. }
    \item{\code{Classification}}{Classification of colon cancer. }
  }
  The following information on the measured genes is saved in the variables 
  \code{Symbol} and \code{Gene.name} of the \code{featureData} slot.
  \describe{
    \item{\code{UBC}}{Ubiquitin C. }
    \item{\code{UBB}}{Ubiquitin B. }
    \item{\code{SUI1}}{Putative translation initiation factor. }
    \item{\code{NACA}}{Nascent-polypeptide-associated complex alpha polypeptide. }
    \item{\code{FLJ20030}}{Hypothetical protein FLJ20030. }
    \item{\code{CFL1}}{Cofilin 1 (non-muscle). }
    \item{\code{ACTB}}{Actin, beta. }
    \item{\code{CLTC}}{Clathrin, heavy polypeptide (Hc). }
    \item{\code{RPS13}}{Ribosomal protein S13. }
    \item{\code{RPS23}}{Ribosomal protein S23. }
    \item{\code{GAPD}}{Glyceraldehyde-3-phosphate dehydrogenase. }
    \item{\code{TPT1}}{Tumor protein, translationally controlled 1. }
    \item{\code{TUBA6}}{Tubulin alpha 6. }
  }
  For a detailed annotation see Table 1 in Anderson et al. (2004).
}
\details{
  The genes included in this data set were selected by screening 161 colon 
  sample expression profiles.
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
  data(Colon)
  Colon
  head(exprs(Colon))
  pData(Colon)
  fData(Colon)
}
\keyword{datasets}
