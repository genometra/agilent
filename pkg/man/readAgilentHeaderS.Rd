\name{readAgilentHeaderS}
\alias{readAgilentHeaderS}
\title{Read Agilent headers.}
\usage{
readAgilentHeaderS(files = dir(pattern = ".txt"), samplenames = sub(".txt",
  "", files), verbose = TRUE, includeFileNames = FALSE)
}
\arguments{
  \item{files}{names of the files which header is to be
  read.}

  \item{samplenames}{names to be given to the arrays.}

  \item{verbose}{verbose}

  \item{includeFileNames}{should the filename be included
  into datasets?}
}
\value{
A list of 3 components containing the FEPARAMS, STATS, and
FEATURE TYPE FROM the header of the Agilent files.
}
\description{
Reads Agilent header from one or several txt exported
files.
}
\details{
FEPARAMS and STATS from header of all arrays are stored in
a data.frames. FEATURE TYPE is stored in a character vector
(named fieldType) when it is the same for all samples. If
it is not the same across arrays (or if includeFileNames =
TRUE) a data.frame is used an a warning message issued.
Agilent \strong{field types} are coded in \R standards.
}
\examples{
setwd (file.path (system.file ("exampledata", package = "agilent")))
dir ()

ra <- readAgilentHeaderS ()
names (ra)
dim (ra$feparams)
ra$feparams[,1:3]
class (ra[[3]])
ra[[3]]

ra <- readAgilentHeaderS (files = c("agilent1.txt", "agilent2.txt"), samplenames = c ("first", "second"), includeFileNames = FALSE)
names (ra)
dim (ra$feparams)
ra$feparams[,1:3]
class (ra[[3]])
}
\author{
David Montaner \email{dmontaner@cipf.es}
}
\seealso{
\code{\link{readAgilentHeader}} \code{\link{readAgilent}}
\code{\link{readGPR}}
}
\keyword{agilent}
\keyword{header}
\keyword{read}

