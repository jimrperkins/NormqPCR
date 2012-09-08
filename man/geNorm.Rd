\name{geNorm}
\alias{geNorm}
\alias{geNorm.qPCRBatch}
\docType{data}
\title{ Data set of Vandesompele et al (2002) }
\description{
  This data set was used in Vandesompele et al (2002) to demonstrate normalization 
  of real-time quantitative RT-PCR data by geometric averaging of housekeeping genes.
}
\usage{data(geNorm)}
\format{
  A qPCRBatch object which contains an expression matrix with 85 observations on the following 10 variables which stand for 
  expression data of ten potential reference/housekeeping genes
  \describe{
    \item{\code{ACTB}}{actin, beta}
    \item{\code{B2M}}{beta-2-microglobulin}
    \item{\code{GAPD}}{glyceraldehyde-3-phosphate dehydrogenase}
    \item{\code{HMBS}}{hydroxymethylbilane synthase}
    \item{\code{HPRT1}}{hypoxanthine phosphoribosyltransferase 1}
    \item{\code{RPL13A}}{ribosomal protein L13a}
    \item{\code{SDHA}}{succinate dehydrogenase complex subunit A}
    \item{\code{TBP}}{TATA box binding protein}
    \item{\code{UBC}}{ubiquitin C}
    \item{\code{YWHAZ}}{tyrosine 3-monooxygenase/tryptophan 5-monooxygenase activation protein, zeta polypeptide}
  }
}
\details{
  The row names of this data set indicate the various human tissues which were
  investigated. 
  \describe{
    \item{BM}{9 normal bone-marrow samples}
    \item{POOL}{9 normal human tissues from pooled organs (heart, brain, fetal brain, lung, 
                trachea, kidney, mammary gland, small intestine and uterus)}
    \item{FIB}{20 short-term cultured normal fibroblast samples from different individuals}
    \item{LEU}{13 normal leukocyte samples}
    \item{NB}{34 neuroblastoma cell lines (independently prepared in different labs from different patients)}
  }
}
\source{
  The data set was obtained from
  \url{http://genomebiology.com/content/supplementary/gb-2002-3-7-research0034-s1.txt}
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
\examples{
  data(geNorm)
  str(exprs(geNorm.qPCRBatch))
  sampleNames(geNorm.qPCRBatch)
}
\keyword{datasets}
