##foreg.r
##2009-10-18 dmontaner@cipf.es
##2013-11-29 dmontaner@cipf.es


##' @name foreg
##' @author David Montaner \email{dmontaner@@cipf.es}
##'
##' @aliases foreground
##' @keywords foreground
##' @seealso \code{\link{backg}} \code{\link{flags}} \code{\link{readAgilentHeader}} \code{\link{readGPR}}
##'
##' @title Extract Foreground form an ExpressionSet.
##' 
##' @description Retrieve foreground data from Agilent One Color ExpressionSet.
##' 
##' @details This function access the F matrix in the assayData of an Agilent One Color ExpressionSet.
##'
##' @param object an Agilent One Color ExpressionSet.
##'
##' @return F matrix from the Agilent ExpressionSet
##' 
##' @export

foreg <- function (object) {
  if (!validObject (object)) stop ("object is not a valid ExpressionSet")
  if (!"F" %in% assayDataElementNames (object)) stop ("the ExpressionSet does not have an F assayData element")
  return (assayDataElement (object, "F"))
}

################################################################################


##' @name backg
##' @author David Montaner \email{dmontaner@@cipf.es}
##'
##' @aliases background
##' @keywords background
##' @seealso \code{\link{foreg}} \code{\link{flags}} \code{\link{readAgilentHeader}} \code{\link{readGPR}}
##'
##' @title Extract Background form an ExpressionSet.
##' 
##' @description Retrieve background data from Agilent One Color ExpressionSet.
##' 
##' @details This function access the B matrix in the assayData of an Agilent One Color ExpressionSet.
##' If notfound2zero is FALSE, when B matrix is not found in the assayData, a NULL value is returned.
##' 
##' @param object an Agilent One Color ExpressionSet.
##' @param notfound2zero if TRUE, B matrix is filled with zero if background is not present.
##' @param verbose verbose
##'
##' @return B matrix from the Agilent ExpressionSet
##' 
##' @export

backg <- function (object, notfound2zero = FALSE, verbose = TRUE) {
  if (!validObject (object)) stop ("object is not a valid ExpressionSet")
  if ("B" %in% assayDataElementNames (object)) {
    out <- assayDataElement (object, "B")
  } else {
    if (verbose) warning ("the ExpressionSet does not have a B assayData element")
    if (notfound2zero) {
      if (verbose) warning ("All elements are set to 0")
      out <-  matrix (0, nrow = nrow (object), ncol = ncol (object))
      dimnames (out) <- dimnames (exprs (object))
    } else {
      out <- NULL
    }
  }
  return (out)
}

################################################################################

##' @name flags
##' @author David Montaner \email{dmontaner@@cipf.es}
##'
##' @keywords flags
##' @seealso \code{\link{foreg}} \code{\link{backg}} \code{\link{readAgilentHeader}} \code{\link{readGPR}}
##'
##' @title Extract flagout form an ExpressionSet.
##' 
##' @description Retrieve flagout information data from Agilent One Color ExpressionSet.
##' 
##' @details This function tries to access the flagout matrix in the assayData of an Agilent One Color ExpressionSet.
##' If notfound2zero is FALSE, when B matrix is not found in the assayData, a NULL value is returned.
##' 
##' @param object an Agilent One Color ExpressionSet.
##' @param notfound2zero If TRUE and flagout matrix is not found in the assayData a zero matrix is returned.
##' @param verbose verbose.
##' 
##' @return B matrix from the Agilent ExpressionSet
##' 
##' @export

flags <- function (object, notfound2zero = TRUE, verbose = FALSE) {
  if (!validObject (object)) stop ("object is not a valid ExpressionSet")
  if ("flagout" %in% assayDataElementNames (object)) {
    out <- assayDataElement (object, "flagout")
  } else {
    if (verbose) warning ("the ExpressionSet does not have a flagout assayData element")
    if (notfound2zero) {
      out <-  matrix (0, nrow = nrow (object), ncol = ncol (object))
      dimnames (out) <- dimnames (exprs (object))
    } else {
      out <- NULL
    }
  }
  return (out)
}
