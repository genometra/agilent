##normalizeAgilent.r
##2009-11-11 dmontaner@cipf.es
##2012-05-18 dmontaner@cipf.es: added the scale.noabs method
##2013-11-29 dmontaner@cipf.es

##' @name normalizeAgilent
##' @author David Montaner \email{dmontaner@@cipf.es}
##'
##' @keywords normalize agilent
##' @seealso \code{\link{readAgilent}} \code{\link{readGPR}} \code{\link{backgroundCorrectAgilent}}
##'
##' @title Normalize Agilent.
##' 
##' @description Normalization methods for Agilent one color microarrays.
##' 
##' @details Available methods are: "quantiles", "scale" and "none".
##' 
##' See \code{\link{backgroundCorrect}} in the limma package.
##' 
##' A column "BiologicalFeature" is expected in the fData if nonbiologicalout is TRUE
##' flagout uses the assayData matrix "flagout" if exists. This is a logical or a 0,1 matrix. See \code{flags}
##' "quantiles" method performs as the function normalize.quantiles in the affy/preprocessCore package
##' "scale" method performs as the function normalizeMedianAbsValues in the limma package.
##' "scale.noabs" performs as the function normalizeMedianValues in the limma package (no absolute value is taken).
##' 
##' @param object an Agilent One Color ExpressionSet.
##' @param method normalization  method
##' @param logit should data be log2 transformed after normalization
##' @param nonbiologicalout should Non Biological Features be discarded when estimating correction parameters
##' @param nonbiological2na should Non Biological Features be set to NA
##' @param flagout should FLAGGED spots be discarded when estimating correction parameters
##' @param flag2na should FLAGGED spots be set to NA
##' @param keep should assayData elements other than "exprs" be kept in the output
##' @param verbose verbose
##' @param \dots further arguments NOT IN USE
##' 
##' @return An ExpressionSet containing normalized intensities in the "exprs" slot.
##' 
##' @examples
##' library (Biobase)
##' setwd (file.path (system.file ("exampledata", package = "agilent")))
##' dir ()
##' 
##' ## ra <- readAgilent ()
##' ## ba <- backgroundCorrectAgilent (ra)
##' 
##' rg <- readGPR ()
##' fData (rg)[1:3,]
##' fData (rg)[,"BiologicalFeature"] <- fData (rg)[,"ControlType"] == "ignore"
##' ba <- backgroundCorrectAgilent (rg)
##' 
##' assayDataElementNames (ba)
##' summary (exprs (ba))
##' 
##' norm <- normalizeAgilent (ba)
##' assayDataElementNames (norm)
##' 
##' norm1 <- normalizeAgilent (object = ba, method = "quantiles", logit = FALSE,
##'                            nonbiologicalout = FALSE, nonbiological2na = FALSE,
##'                            flagout = FALSE, flag2na = FALSE,
##'                            keep = TRUE, verbose = TRUE)
##' 
##' table (exprs (norm1) == preprocessCore::normalize.quantiles (exprs (ba)), exclude = NULL)
##' 
##' norm2 <- normalizeAgilent (object = ba, method = "scale", logit = FALSE,
##'                            nonbiologicalout = FALSE, nonbiological2na = FALSE,
##'                            flagout = FALSE, flag2na = FALSE,
##'                            keep = TRUE, verbose = TRUE)
##' 
##' table (exprs (norm2) == limma::normalizeMedianAbsValues (exprs (ba)), exclude = NULL)
##' 
##' norm3 <- normalizeAgilent (object = ba, method = "quantiles")
##' summary (exprs (norm3))
##' summary (exprs (norm1))
##' 
##' norm4 <- normalizeAgilent (object = ba, method = "scale")
##' summary (exprs (norm4))
##' summary (exprs (norm2))
##' 
##' norm5 <- normalizeAgilent (object = ba, method = "quantiles", logit = FALSE,
##'                            nonbiologicalout = TRUE, nonbiological2na = TRUE,
##'                            flagout = FALSE, flag2na = FALSE,
##'                            keep = TRUE, verbose = TRUE)
##' 
##' table (fData(ba)[,"BiologicalFeature"] == fData(norm5)[,"BiologicalFeature"], exclude = NULL)
##' biological <- fData(ba)[,"BiologicalFeature"]
##' table (biological, exclude = NULL)
##' 
##' summary (exprs (norm5[!biological,]))
##' summary (exprs (norm5[biological,]))
##' summary (preprocessCore::normalize.quantiles (exprs (ba[biological,])))
##' summary (exprs (norm5[biological,]) -
##'          preprocessCore::normalize.quantiles (exprs (ba[biological,]))) #similar
##' 
##' norm6 <- normalizeAgilent (object = ba, method = "scale", logit = FALSE,
##'                            nonbiologicalout = TRUE, nonbiological2na = TRUE,
##'                            flagout = FALSE, flag2na = FALSE,
##'                            keep = TRUE, verbose = TRUE)
##' 
##' table (fData(ba)[,"BiologicalFeature"] == fData(norm6)[,"BiologicalFeature"], exclude = NULL)
##' biological <- fData(ba)[,"BiologicalFeature"]
##' table (biological, exclude = NULL)
##' 
##' summary (exprs (norm6[!biological,]))
##' summary (exprs (norm6[biological,]))
##' summary (limma::normalizeMedianAbsValues (exprs (ba[biological,])))
##' table (exprs (norm6[biological,]) ==
##'        limma::normalizeMedianAbsValues (exprs (ba[biological,]))) #similar
##'
##' @import preprocessCore
## @importFrom preprocessCore normalize.quantiles.determine.target normalize.quantiles.use.target
##' 
##' @export

normalizeAgilent <- function (object, method = "quantiles", logit = FALSE,
                                      nonbiologicalout = TRUE, nonbiological2na = FALSE,
                                      flagout = TRUE, flag2na = FALSE,
                                      keep = FALSE, verbose = TRUE, ...) {#2009-11-11

  ##checking methods
  methods <- c ("quantiles", "scale", "scale.noabs", "none") #"qrobust"
  if (!method %in% methods) stop (paste (method, "is not a valid method"))

  
  ##checking Agilent One Color ExpressionSet object
  if (!validObject (object)) stop (paste (object, "is not a valid ExpressionSet"))
  ##   if (!"F" %in% assayDataElementNames (object)) stop (paste (object, "does not have an F assayData element"))
  ##   if (!"B" %in% assayDataElementNames (object)) warning ("object does not have an B assayData element. Some methods may fail.")
  
  if (flagout | flag2na) {
    if (!"flagout" %in% assayDataElementNames (object)) warning ("object does not have a flagout assayData element.")
  }
  
  ##taking out non "BiologicalFeature"
  if (nonbiologicalout | nonbiological2na) {
    if (!"BiologicalFeature" %in% colnames (fData(object))) {
      stop (paste (object, "does not have a BiologicalFeature column in pData. Nonbiological features cannot be excluded"))
    }
  } 
  
  
  ##NORMALIZING ----------------------------------------------------------------

  if (verbose) message (paste ("\nnormalization  method:",  method))
  
  pred <- train <- exprs (object)
  ##train
  if (flagout | nonbiologicalout) {
    excludedByFlag <- flagout & flags (object) #nos da la matriz de la dimension necesaria
    excludedNonBiological <- FALSE
    if (nonbiologicalout) {
      excludedNonBiological <- nonbiologicalout & !fData(object)[,"BiologicalFeature"]
  }
    exclude <- excludedByFlag | excludedNonBiological
    train[exclude] <- NA
  }
  ##pred
  if (flag2na | nonbiological2na) {
    excludedByFlag <- flag2na & flags (object) #nos da la matriz de la dimension necesaria
    excludedNonBiological <- FALSE
    if (nonbiological2na) {
      excludedNonBiological <- nonbiological2na & !fData(object)[,"BiologicalFeature"]
    }
    exclude <- excludedByFlag | excludedNonBiological
    pred[exclude] <- NA
  }
  
  if (method == "quantiles") {#similar to normalize.quantiles from the affy/preprocessCore package
    ##require (preprocessCore)  
    tgts <- normalize.quantiles.determine.target (train)
    exprs (object) <- normalize.quantiles.use.target (pred, tgts)
    rownames (exprs (object)) <- rownames (pred)
    colnames (exprs (object)) <- colnames (pred)
  }
  
  if (method == "scale") {#as in normalizeMedianAbsValues from the limma package
    ##require (limma)  
    if (ncol (object) > 1) { #with just one column scale does nothing
      cmed <- log (apply (abs (train), 2, median, na.rm = TRUE))  ##as in normalizeMedianAbsValues
      #cmed <- log (apply (train, 2, median, na.rm = TRUE))       ##as in normalizeMedianValues
      cmed <- exp (cmed - mean(cmed))
      exprs (object) <- t ( t (pred) / cmed)
    }
  }

  if (method == "scale.noabs") {#as in normalizeMedianAbsValues from the limma package
    ##require (limma)  
    if (ncol (object) > 1) { #with just one column scale does nothing
      #cmed <- log (apply (abs (train), 2, median, na.rm = TRUE))  ##as in normalizeMedianAbsValues
      cmed <- log (apply (train, 2, median, na.rm = TRUE))         ##as in normalizeMedianValues
      cmed <- exp (cmed - mean(cmed))
      exprs (object) <- t ( t (pred) / cmed)
    }
  }



  
  ##log ------------------------------------------------------------------------

  if (logit) {
    exprs (object) <- log2 (exprs (object))
    if (verbose) message ("\nlog2 transformation done after normalization")
  }

  ##keep: eliminate some assayDataElements -------------------------------------
  if (!keep) {
    elementos <- setdiff (assayDataElementNames (object), "exprs")
    if (verbose) message (paste ( c("\nElements: {", elementos, "} will be removed from the assayData"), collapse = " "))
    for (e in elementos) {
      object <- assayDataElementReplace (object, e, NULL)
    }
  }

  ##chequeamos que el objeto de clase "ExpressionSet" es valido ----------------
  validObject (object)
  
  ##SALIDA
  return (object)
}


####################################################################################################################
####################################################################################################################


## library (agilent,  lib.loc = "/home/david/programas/mislibrerias/agilent/agilent_instalacion_local"); packageDescription ("agilent", fields = "Version") #

## setwd (file.path (system.file ("exampledata", package = "agilent")))
## dir ()

## ## ra <- readAgilent ()
## ## ba <- backgroundCorrectAgilent (ra)

## rg <- readGPR ()
## fData (rg)[1:3,]
## fData (rg)[,"BiologicalFeature"] <- fData (rg)[,"ControlType"] == "ignore"
## ba <- backgroundCorrectAgilent (rg)

## assayDataElementNames (ba)
## summary (exprs (ba))

## norm <- normalizeAgilent (ba)
## assayDataElementNames (norm)

## norm1 <- normalizeAgilent (object = ba, method = "quantiles", logit = FALSE,
##                            nonbiologicalout = FALSE, nonbiological2na = FALSE,
##                            flagout = FALSE, flag2na = FALSE,
##                            keep = TRUE, verbose = TRUE)

## table (exprs (norm1) == normalize.quantiles (exprs (ba)), exclude = NULL)

## norm2 <- normalizeAgilent (object = ba, method = "scale", logit = FALSE,
##                            nonbiologicalout = FALSE, nonbiological2na = FALSE,
##                            flagout = FALSE, flag2na = FALSE,
##                            keep = TRUE, verbose = TRUE)

## table (exprs (norm2) == normalizeMedianAbsValues (exprs (ba)), exclude = NULL)

## ###

## norm3 <- normalizeAgilent (object = ba, method = "quantiles")
## summary (exprs (norm3))
## summary (exprs (norm1))

## norm4 <- normalizeAgilent (object = ba, method = "scale")
## summary (exprs (norm4))
## summary (exprs (norm2))


## norm5 <- normalizeAgilent (object = ba, method = "quantiles", logit = FALSE,
##                            nonbiologicalout = TRUE, nonbiological2na = TRUE,
##                            flagout = FALSE, flag2na = FALSE,
##                            keep = TRUE, verbose = TRUE)

## table (fData(ba)[,"BiologicalFeature"] == fData(norm5)[,"BiologicalFeature"], exclude = NULL)
## biological <- fData(ba)[,"BiologicalFeature"]
## table (biological, exclude = NULL)

## summary (exprs (norm5[!biological,]))
## summary (exprs (norm5[biological,]))
## summary (normalize.quantiles (exprs (ba[biological,])))
## summary (exprs (norm5[biological,]) - normalize.quantiles (exprs (ba[biological,]))) #similar


## norm6 <- normalizeAgilent (object = ba, method = "scale", logit = FALSE,
##                            nonbiologicalout = TRUE, nonbiological2na = TRUE,
##                            flagout = FALSE, flag2na = FALSE,
##                            keep = TRUE, verbose = TRUE)

## table (fData(ba)[,"BiologicalFeature"] == fData(norm6)[,"BiologicalFeature"], exclude = NULL)
## biological <- fData(ba)[,"BiologicalFeature"]
## table (biological, exclude = NULL)

## summary (exprs (norm6[!biological,]))
## summary (exprs (norm6[biological,]))
## summary (normalizeMedianAbsValues (exprs (ba[biological,])))
## table (exprs (norm6[biological,]) == normalizeMedianAbsValues (exprs (ba[biological,]))) #similar
