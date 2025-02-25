% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deriveAgilentFeatures.r
\name{deriveAgilentFeatures}
\alias{deriveAgilentFeatures}
\title{Derive Agilent Feature Data.}
\usage{
deriveAgilentFeatures(data, log.relative.conc = NULL, verbose = TRUE)
}
\arguments{
\item{data}{a data.frame of feature Data obtained from readAgilent}

\item{log.relative.conc}{named vector of SpikeIns log relative concentration.}

\item{verbose}{verbose}
}
\value{
A data.frame with derived features.
}
\description{
Describes Agilent Feature data.
}
\details{
Derives some feature data for the spots of the array.
Agilent features "SubTypeName" and "ControlType" are combined into a new "FeatureType".
Notice that Agilent features "SubTypeName" and "SubTypeMask" are equivalent.
Only spots tagged as "BiologicalFeature" in the "FeatureType" are suitable for analyzing.
Spots tagged as "E1A" in the "FeatureType" are SpikeIns.
log.relative.conc is a named vector describing
1-color Agilent SpikeIns Signal Statistics Log (Relative Conc.)
Concentrations form Agilent Design ID 014850 are assumed by default.
}
\seealso{
\code{\link{readAgilent}} \code{\link{readGPR}}
}
\author{
David Montaner \email{dmontaner@cipf.es}
}
\keyword{agilent}
\keyword{features}
