\name{stabMeasureRho}
\alias{stabMeasureRho}
\alias{stabMeasureRho-methods}
\alias{stabMeasureRho,qPCRBatch-method}
\alias{stabMeasureRho,matrix-method}
\alias{stabMeasureRho,x-method}
\title{ Gene expression stability value rho }
\description{
  Computation of the gene expression stability value rho for real-time 
  quantitativ RT-PCR data. For more details we refer to 
  Andersen et al. (2004). 
}
\usage{
  stabMeasureRho(x,\dots)

  \S4method{stabMeasureRho}{x}(x, group, log = TRUE, na.rm = TRUE, returnAll = FALSE)
}
\arguments{
  \item{x}{ matrix containing real-time quantitative 
  		  RT-PCR data, or qPCRBatch object }
  \item{\dots}{ Extra arguments, detailed below }
  \item{group}{ grouping factor, either a factor vector or a phenoData column called "Group" }  
  \item{log}{ logical: is data on log-scale }
  \item{na.rm}{ a logical value indicating whether \code{NA} values should be
          stripped before the computation proceeds. }
  \item{returnAll}{ logical, return additional information. }
}
\details{
  The gene expression stability value rho is computed. For more details see 
  Andersen et al. (2004).
}
\value{
  numeric vector with gene expression stability values
  
  If \code{returnAll == TRUE} a list with the following components is returned
  \item{rho }{ stability measure rho of Andersen et al. (2004) }
  \item{d }{ used by \code{selectHKs} }
  \item{v }{ used by \code{selectHKs} }
}
\references{ 
  Claus Lindbjerg Andersen, Jens Ledet Jensen and Torben Falck Orntoft (2004).
  Normalization of Real-Time Quantitative Reverse Transcription-PCR Data: 
  A Model-Based Variance Estimation Approach to Identify Genes Suited for
  Normalization, Applied to Bladder and Colon Cancer Data Sets. 
  CANCER RESEARCH 64, 5245-5250, August 1, 2004.
  \url{http://cancerres.aacrjournals.org/cgi/content/full/64/15/5245}
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{selectHKs}}
\examples{
  data(Colon)
  group <- pData(Colon.qPCRBatch)[,"Group"]
  res.Colon <- stabMeasureRho(Colon.qPCRBatch, group = group,
                            log = FALSE)
  sort(res.Colon) # cf. Table 3 in Andersen et al (2004)
}
\keyword{data}
