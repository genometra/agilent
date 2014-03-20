\name{readAgilent}
\alias{readAgilent}
\title{Read Agilent.}
\usage{
  readAgilent(files = character(0), samplenames = NULL,
    sampleinfo = NULL, foreground.column = "gMeanSignal",
    background.column = "gBGMedianSignal",
    other.columns = c("gProcessedSignal", "IsManualFlag"),
    feature.columns = character(0), derivefeatures = TRUE,
    verbose = TRUE, ...)
}
\arguments{
  \item{files}{names of the files to be read.}

  \item{samplenames}{names to be given to the arrays.}

  \item{sampleinfo}{data frame with columns 'filename',
  'samplename' and any other sample information.}

  \item{foreground.column}{column name of the foreground
  intensities.}

  \item{background.column}{column name of the background
  intensities.}

  \item{other.columns}{some other columns to be read.}

  \item{feature.columns}{some other columns to be read as
  featureData.}

  \item{derivefeatures}{should additional features be
  computed.}

  \item{verbose}{verbose}

  \item{\dots}{further arguments to be passed to
  deriveAgilentFeatures.}
}
\value{
  An ExpressionSet containing foreground (F) and background
  (B) matrices within the assayData
}
\description{
  Reads Agilent raw data txt exported files.
}
\details{
  exprs in assayData is set to NA

  sampleinfo is used rather than files or samplenames

  sampleinfo must have columns named 'filename' and
  'samplename'

  If background.columncolumn = NULL no B matrix is set in
  assayData. If all feature.columns are the same a single
  column is passed to the featureData. If not all of them
  are equal, the information stays in the assayData

  Default feature.columns are: "ProbeName", "GeneName",
  "SystematicName", "FeatureNum", "Row", "Col",
  "SubTypeMask", "SubTypeName", "ProbeUID", "ControlType",
  "Description"

  Extra columns read by default are: "gProcessedSignal" and
  "IsManualFlag". See the \code{other.columns} parameter.
}
\examples{
setwd (file.path (system.file ("exampledata", package = "agilent")))
dir ()

ra <- readAgilent ()
pData (ra)
fData (ra)[1:5,]
assayDataElementNames (ra)

sinfo <- cbind (filename = c("agilent2.txt", "agilent1.txt"),
                samplename = c("first", "second"), more = 1:2)
ra <- readAgilent (sampleinfo = sinfo)
pData (ra)

ra <- readAgilent (background.column = NULL, other.columns = NULL, feature.columns = NULL, verbose = FALSE)
pData (ra)
fData (ra)[1:5,]
assayDataElementNames (ra)
}
\author{
  David Montaner \email{dmontaner@cipf.es}
}
\seealso{
  \code{\link{readAgilentHeader}} \code{\link{readGPR}}
}
\keyword{agilent}
\keyword{read}

