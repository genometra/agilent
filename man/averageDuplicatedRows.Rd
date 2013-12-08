\name{averageDuplicatedRows}
\alias{averageDuplicatedRows}
\alias{insert}
\title{Average duplicated rows}
\usage{
  averageDuplicatedRows(mat, ids, original.rows = FALSE,
    summarystat = "mean", na.rm = TRUE, verbose = TRUE)
}
\arguments{
  \item{mat}{a numerical matrix.}

  \item{ids}{character vector of length equal to the number
  of rows of \code{mat}.}

  \item{original.rows}{logical indicating whether the
  mapping between \code{ids} and original row position is
  to be kept.}

  \item{summarystat}{the summary statistic: "mean" or
  "median"}

  \item{na.rm}{should NA values be removed before
  computation?}

  \item{verbose}{verbose}
}
\value{
  a matrix of averaged rows. \code{ids} are used as row
  names.
}
\description{
  Averages duplicated rows in a matrix.
}
\details{
  The function averages rows in a matrix that are indexed
  by the same ids.

  If original.rows is \code{TRUE} a length 2 list is
  returned. The first element contains the row averaged
  values. The second one keeps the link between ids and
  original row position.
}
\examples{
X <- data.frame(go=c("a","b","a","b"), exp1=c(3,4,5,7),
                exp2=c(8,9,15,20), stringsAsFactors = FALSE)
X
M <- as.matrix (X[,c("exp1", "exp2")])
averageDuplicatedRows (mat = M, id = X[,"go"])
averageDuplicatedRows (mat = M, id = X[,"go"], original.rows = TRUE)
}
\author{
  David Montaner \email{dmontaner@cipf.es}
}
\keyword{insert}
\keyword{tiddler}

