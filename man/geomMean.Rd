\name{geomMean}
\alias{geomMean}
\title{ Geometric Mean }
\description{
  Computation of the geometric mean.
}
\usage{
geomMean(x, na.rm = TRUE)
}
\arguments{
  \item{x}{ numeric vector of non-negative Reals }
  \item{na.rm}{ a logical value indicating whether \code{NA} values should be
          stripped before the computation proceeds.  }
}
\details{
  The computation of the geometric mean is done via \code{prod(x)^(1/length(x))}.
}
\value{
  geometric mean
}
\references{ 
  Perkins, JR, Dawes, JM, McMahon, SB, Bennett, DL, Orengo, C, Kohl, M (2012). 
  ReadqPCR and NormqPCR: R packages for the reading, quality checking and 
  normalisation of RT-qPCR quantification cycle (Cq) data. 
  \emph{BMC Genomics}, \bold{13}, 1:296.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\note{
A first version of this function appeared in package SLqPCR.
}
%\seealso{}
\examples{
x <- rlnorm(100)
geomMean(x)
}
\keyword{data}
