\name{Bladder}
\alias{Bladder}
\alias{Bladder.qPCRBatch}
\docType{data}
\title{ Data set of Andersen et al (2004) }
\description{
  This data set was used in Andersen et al (2004) to demonstrate normalization 
  of real-time quantitative RT-PCR data by geometric averaging of housekeeping genes.
}
\usage{data(Bladder)}
\format{
  A qPCRBatch object which contains an expression matrix with 28 observations on the following 14 variables which stand for 
  expression data potential reference/housekeeping genes
  \describe{
    \item{\code{Grade}}{Grade of bladder cancer. }
    \item{\code{Sample.no.}}{sample number. }
    \item{\code{ATP5B}}{ATP synthase, H+ transporting, mitochondrial F1 complex, beta polypetide. }
    \item{\code{HSPCB}}{Heat shock 90-kDa protein 1, beta. }
    \item{\code{S100A6}}{S100 calcium-binding protein A6 (calcylin). }
    \item{\code{FLOT2}}{Flotillin 2. }
    \item{\code{TEGT}}{Testis enhanced gene transcript (BAX inhibitor 1). }
    \item{\code{UBB}}{Ubiquitin B. }
    \item{\code{TPT1}}{Tumor protein, translationally controlled 1. }
    \item{\code{CFL1}}{Cofilin 1 (non-muscle). }
    \item{\code{ACTB}}{Actin, beta. }
    \item{\code{RPS13}}{Ribosomal protein S13. }
    \item{\code{RPS23}}{Ribosomal protein S23. }
    \item{\code{GAPD}}{Glyceraldehyde-3-phosphate dehydrogenase. }
    \item{\code{UBC}}{Ubiquitin C. }
    \item{\code{FLJ20030}}{Hypothetical protein FLJ20030. }
  }
  It also contains a "group" column in the pData slot, which gives information on the different sample classes,
  necessary for the NormFinder algorithm
}
\details{
  The genes included in this data set were selected by screening 99 bladder 
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
  data(Bladder)
  head(exprs(Bladder.qPCRBatch))
  pData(Bladder.qPCRBatch)
}
\keyword{datasets}
