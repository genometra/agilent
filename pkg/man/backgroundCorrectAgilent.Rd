\name{backgroundCorrectAgilent}
\alias{backgroundCorrectAgilent}
\title{Background Correct Agilent.}
\usage{
backgroundCorrectAgilent(object, method = "rma", logit = TRUE,
  nonbiologicalout = TRUE, nonbiological2na = FALSE, flagout = TRUE,
  flag2na = FALSE, keep = TRUE, verbose = TRUE, ...)
}
\arguments{
  \item{object}{an Agilent One Color ExpressionSet.}

  \item{method}{background correction method}

  \item{logit}{should data be log2 transformed after
  background correction}

  \item{nonbiologicalout}{should Non Biological Features be
  discarded when estimating correction parameters}

  \item{nonbiological2na}{should Non Biological Features be
  set to NA}

  \item{flagout}{should FLAGGED spots be discarded when
  estimating correction parameters}

  \item{flag2na}{should FLAGGED spots be set to NA}

  \item{keep}{should assayData elements other than "exprs"
  be kept in the output}

  \item{verbose}{verbose}

  \item{\dots}{further arguments NOT IN USE}
}
\value{
An ExpressionSet containing background corrected features
in the "exprs" slot.
}
\description{
Background correction methods for Agilent one color
microarrays.
}
\details{
Available methods are: "agilent", "none", "subtract",
"half", "normexp", "rma", "normexpnobg", "rmanobg"
"normexp" and "rma" are applied after background
subtraction as in limma implementation. "normexpnobg" and
"rmanobg" are applied without background subtraction. This
may be suitable for miRNA array normalization as described
in Lopez-Romero et al. (2010)
\url{http://www.biomedcentral.com/1756-0500/3/18}

The method "agilent" returns the column "gProcessedSignal"
read from the raw data file. No extra processing is done
with this column. See \code{\link{readAgilent}} for details
on how to read this column.

See \code{\link{backgroundCorrect}} in the limma package.

nonbiologicalout and flagout have an effect only when using
methods "normexp" or "rma.

A column "BiologicalFeature" is expected in the fData if
nonbiologicalout is TRUE flagout uses the assayData matrix
"flagout" if exists. This is a logical or a 0,1 matrix. See
\code{flags}

Depends on the function \code{normexp.fit} in the library
\code{limma}. the function \code{normexp.fit} performs a
\code{require(affy)} and thus the library \code{affy} gets
attached to the path.
}
\examples{
setwd (file.path (system.file ("exampledata", package = "agilent")))
dir ()

ra <- readAgilent ()
assayDataElementNames (ra)

ba <- backgroundCorrectAgilent (ra)
ba

rg <- readGPR ()
assayDataElementNames (rg)
fData (rg)[1:3,]
fData (rg)[,"BiologicalFeature"] <- fData (rg)[,"ControlType"] == "ignore"

bg <- backgroundCorrectAgilent (rg)
summary (exprs (bg))

bg <- backgroundCorrectAgilent (rg, method = "normexp")
summary (exprs (bg))

bg <- backgroundCorrectAgilent (rg, nonbiological2na = TRUE)
summary (exprs (bg))
}
\author{
David Montaner \email{dmontaner@cipf.es}
}
\seealso{
\code{\link{readAgilent}} \code{\link{readGPR}}
}
\keyword{agilent}
\keyword{background}
\keyword{correction}

