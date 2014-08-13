\name{NormqPCR-package}
\alias{NormqPCR-package}
\alias{NormqPCR}
\docType{package}
\title{
  Functions for normalisation of RT-qPCR data.
}
\description{
Functions for normalisation of real-time quantitative PCR data.
}
\details{
\tabular{ll}{
Package: \tab NormqPCR\cr
Type: \tab Package\cr
Version: \tab 1.7.1\cr
Date: \tab 2014-08-13\cr
Depends: \tab R(>= 2.14.0), stats, RColorBrewer, Biobase, methods, ReadqPCR, qpcR \cr
Imports: \tab ReadqPCR \cr
biocViews: \tab MicrotitrePlateAssay, GeneExpression, qPCR \cr
License: \tab LGPL-3\cr
LazyLoad: \tab yes \cr
LazyData: \tab yes \cr

}
require(NormqPCR)
}
\author{ Matthias Kohl, James Perkins, Nor Izayu Abdul Rahman

Maintainer: James Perkins \email{j.perkins@ucl.ac.uk}}
\references{
  Perkins, JR, Dawes, JM, McMahon, SB, Bennett, DL, Orengo, C, Kohl, M (2012). 
  ReadqPCR and NormqPCR: R packages for the reading, quality checking and 
  normalisation of RT-qPCR quantification cycle (Cq) data. 
  \emph{BMC Genomics}, \bold{13}, 1:296.

  Jo Vandesompele, Katleen De Preter, Filip Pattyn et al. (2002). Accurate 
  normalization of real-time quantitative RT-PCR data by geometric averaging 
  of multiple internal control genes. 
  Genome Biology 2002. 3(7):research0034.1-0034.11.
  \url{http://genomebiology.com/2002/3/7/research/0034/}
  
  Claus Lindbjerg Andersen, Jens Ledet Jensen and Torben Falck Orntoft (2004).
  Normalization of Real-Time Quantitative Reverse Transcription-PCR Data: 
  A Model-Based Variance Estimation Approach to Identify Genes Suited for
  Normalization, Applied to Bladder and Colon Cancer Data Sets. 
  CANCER RESEARCH 64, 5245-5250, August 1, 2004.
  \url{http://cancerres.aacrjournals.org/cgi/content/full/64/15/5245}

  Kenneth Livak, Thomase Schmittgen (2001).
  Analysis of Relative Gene Expression 
  Data Using Real-Time Quantitative PCR and the 2^ddCt Method.
  Methods 25, 402-408, 2001
  \url{http://www.ncbi.nlm.nih.gov/pubmed/11846609}
}
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
## some examples are given in the vignette
\dontrun{
library(NormqPCR)
vignette("NormqPCR")
}
}
\keyword{package}
