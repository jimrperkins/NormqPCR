\name{selectHKs}
\alias{selectHKs}
\alias{selectHKs-methods}
\alias{selectHKs,qPCRBatch-method}
\title{ Selection of reference/housekeeping genes }
\description{
  This function can be used to determine a set of reference/housekeeping (HK) 
  genes for gene expression experiments
}
\usage{
selectHKs(qPCRBatch, \dots)

\S4method{selectHKs}{qPCRBatch}(qPCRBatch, group, method = "geNorm", minNrHKs = 2, log = TRUE, Symbols, 
	      trace = TRUE, na.rm = TRUE)
}
\arguments{
  \item{qPCRBatch}{ qPCRBatch, containing the data (expression matrix) in the exprs slot }
  \item{\dots}{ Extra arguments, detailed below }
  \item{group}{ optional factor not used by all methods, hence may be missing }
  \item{method}{ method to compute most stable genes }
  \item{minNrHKs}{ minimum number of HK genes that should be considered }
  \item{log}{ logical: is data on log-scale }
  \item{Symbols}{ gene symbols }
  \item{trace}{ logical, print additional information }
  \item{na.rm}{ a logical value indicating whether \code{NA} values should be
          stripped before the computation proceeds.  }
}
\details{
  This function can be used to determine a set of reference/housekeeping (HK) genes
  for gene expression experiments. The default method \code{"geNorm"}
  was proposed by Vandesompele et al. (2002).

  Currently, the geNorm method by Vandesompele et al. (2002) and the NormFinder
  method of Andersen et al. (2004) are implemented.

  Vandesompele et al. (2002) propose a cut-off value of 0.15 for the 
  pairwise variation. Below this value the inclusion of an additional
  housekeeping gene is not required. 
}
\value{
  If \code{method = "geNorm"} a list with the following components is
  returned
  \item{ranking }{ ranking of genes from best to worst where the two most
  		   stable genes cannot be ranked }
  \item{variation }{ pairwise variation during stepwise selection }
  \item{meanM }{ average expression stability M }

  If \code{method = "NormFinder"} a list with the following components is
  returned
  \item{ranking }{ ranking of genes from best to worst where the two most
         stable genes cannot be ranked }
  \item{rho }{ stability measure rho of Andersen et al. (2004) }
}
\references{ 
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
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{}
\examples{
  data(geNorm)
  tissue <- as.factor(c(rep("BM", 9), rep("FIB", 20), rep("LEU", 13),
    rep("NB", 34), rep("POOL", 9)))

    res.BM <- selectHKs(geNorm.qPCRBatch[,tissue == "BM"], method = "geNorm", 
  Symbols = featureNames(geNorm.qPCRBatch), minNrHK = 2, log = FALSE)
}
\keyword{data}
