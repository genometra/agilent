\name{backg}
\alias{backg}
\alias{background}
\title{Extract Background form an ExpressionSet.}
\usage{
  backg(object, notfound2zero = FALSE, verbose = TRUE)
}
\arguments{
  \item{object}{an Agilent One Color ExpressionSet.}

  \item{notfound2zero}{if TRUE, B matrix is filled with
  zero if background is not present.}

  \item{verbose}{verbose}
}
\value{
  B matrix from the Agilent ExpressionSet
}
\description{
  Retrieve background data from Agilent One Color
  ExpressionSet.
}
\details{
  This function access the B matrix in the assayData of an
  Agilent One Color ExpressionSet. If notfound2zero is
  FALSE, when B matrix is not found in the assayData, a
  NULL value is returned.
}
\author{
  David Montaner \email{dmontaner@cipf.es}
}
\seealso{
  \code{\link{foreg}} \code{\link{flags}}
  \code{\link{readAgilentHeader}} \code{\link{readGPR}}
}
\keyword{background}

