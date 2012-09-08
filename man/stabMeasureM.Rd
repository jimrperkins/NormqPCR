\name{stabMeasureM}
\alias{stabMeasureM}
\title{ Gene expression stability value M }
\description{
  Computation of the gene expression stability value M for real-time 
  quantitativ RT-PCR data. For more details we refer to 
  Vandesompele et al. (2002). 
}
\usage{
stabMeasureM(x, log = TRUE, na.rm = TRUE)
}
\arguments{
  \item{x}{ matrix or data.frame containing real-time quantitative 
  		  RT-PCR data }
  \item{log}{ logical: is data on log-scale }
  \item{na.rm}{ a logical value indicating whether \code{NA} values should be
          stripped before the computation proceeds. }
}
\details{
  The gene expression stability value M is defined as the average pairwise
  normalization factor; i.e., one needs to specify data from at least two
  genes. For more details see Vandesompele et al. (2002). Note this 
  dispatches on a transposed expression matrix, not a qPCRBatch
  object since it is only called from within the selectHKs method.

}
\value{
  numeric vector with gene expression stability values 
}
\references{ 
  Jo Vandesompele, Katleen De Preter, Filip Pattyn et al. (2002). Accurate 
  normalization of real-time quantitative RT-PCR data by geometric averaging 
  of multiple internal control genes. 
  Genome Biology 2002. 3(7):research0034.1-0034.11.
  \url{http://genomebiology.com/2002/3/7/research/0034/}

  Perkins, JR, Dawes, JM, McMahon, SB, Bennett, DL, Orengo, C, Kohl, M (2012). 
  ReadqPCR and NormqPCR: R packages for the reading, quality checking and 
  normalisation of RT-qPCR quantification cycle (Cq) data. 
  \emph{BMC Genomics}, \bold{13}, 1:296.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{selectHKs}}
\examples{
  data(geNorm)
  tissue <- as.factor(c(rep("BM", 9),  rep("FIB", 20), rep("LEU", 13),
                    rep("NB", 34), rep("POOL", 9)))
  res.BM <- selectHKs(geNorm.qPCRBatch[,tissue == "BM"], method = "geNorm", 
                    Symbols = featureNames(geNorm.qPCRBatch), minNrHK = 2, log = FALSE)
}
\keyword{data}
