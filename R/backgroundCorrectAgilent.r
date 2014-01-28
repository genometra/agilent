##backgroundCorrectAgilent.r
##2009-10-20 dmontaner@cipf.es
##2013-11-29 dmontaner@cipf.es


##' @name backgroundCorrectAgilent
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @keywords agilent background correction
##' @seealso \code{\link{readAgilent}} \code{\link{readGPR}}
##'
##' @title Background Correct Agilent.
##' 
##' @description Background correction methods for Agilent one color microarrays.
##'
##' @details Available methods are: "agilent", "none", "subtract", "half", "normexp", "rma", "normexpnobg", "rmanobg"
##' "normexp" and "rma" are applied after background subtraction as in limma implementation.
##' "normexpnobg" and "rmanobg" are applied without background subtraction.
##' This may be suitable for miRNA array normalization as described in Lopez-Romero et al. (2010)
##' \url{http://www.biomedcentral.com/1756-0500/3/18}
##'
##' The method "agilent" returns the column "gProcessedSignal" read from the raw data file.
##' No extra processing is done with this column.
##' See \code{\link{readAgilent}} for details on how to read this column.
##' 
##' See \code{\link{backgroundCorrect}} in the limma package.
##' 
##' nonbiologicalout and flagout have an effect only when using methods "normexp" or "rma.
##'
##' A column "BiologicalFeature" is expected in the fData if nonbiologicalout is TRUE
##' flagout uses the assayData matrix "flagout" if exists. This is a logical or a 0,1 matrix. See \code{flags}
##'
##' Depends on the function \code{normexp.fit} in the library \code{limma}.
##' the function \code{normexp.fit} performs a \code{require(affy)}
##' and thus the library \code{affy} gets attached to the path.
##' 
##' @param object an Agilent One Color ExpressionSet.
##' @param method background correction method
##' @param logit should data be log2 transformed after background correction
##' @param nonbiologicalout should Non Biological Features be discarded when estimating correction parameters
##' @param nonbiological2na should Non Biological Features be set to NA
##' @param flagout should FLAGGED spots be discarded when estimating correction parameters
##' @param flag2na should FLAGGED spots be set to NA
##' @param keep should assayData elements other than "exprs" be kept in the output
##' @param verbose verbose
##' @param \dots further arguments NOT IN USE
##' 
##' @return An ExpressionSet containing background corrected features in the "exprs" slot.
##' 
##' @examples
##' setwd (file.path (system.file ("exampledata", package = "agilent")))
##' dir ()
##'
##' ra <- readAgilent ()
##' assayDataElementNames (ra)
##'
##' ba <- backgroundCorrectAgilent (ra)
##' ba
##'
##' rg <- readGPR ()
##' assayDataElementNames (rg)
##' fData (rg)[1:3,]
##' fData (rg)[,"BiologicalFeature"] <- fData (rg)[,"ControlType"] == "ignore"
##'
##' bg <- backgroundCorrectAgilent (rg)
##' summary (exprs (bg))
##' 
##' bg <- backgroundCorrectAgilent (rg, method = "normexp")
##' summary (exprs (bg))
##' 
##' bg <- backgroundCorrectAgilent (rg, nonbiological2na = TRUE)
##' summary (exprs (bg))
##'
##  @import limma ##WHY IS THIS NOT WORKING
##' @importFrom limma normexp.fit normexp.signal
##' 
##' @export

backgroundCorrectAgilent <- function (object, method = "rma", logit = TRUE,
                                      nonbiologicalout = TRUE, nonbiological2na = FALSE,
                                      flagout = TRUE, flag2na = FALSE,
                                      keep = TRUE, verbose = TRUE, ...) {#2010-10-08

  ##checking methods
  ##methods <- c ("agilent", "none", "subtract", "half", "normexp", "rma")
  methods <- c ("agilent", "none", "subtract", "half", "normexp", "rma", "normexpnobg", "rmanobg")
  if (!method %in% methods) stop (paste (method, "is not a valid method"))

  
  ##checking Agilent One Color ExpressionSet object
  if (!validObject (object)) stop ("object is not a valid ExpressionSet")
  if (!"F" %in% assayDataElementNames (object)) stop ("object does not have an F assayData element")
  if (!"B" %in% assayDataElementNames (object)) warning ("object does not have an B assayData element. Some methods may fail.")

  if (method == "agilent") {
    if (!"gProcessedSignal" %in% assayDataElementNames (object))
      stop ("object does not have a gProcessedSignal assayData element. Background correction method agilent cannot be used.")
  }
  
  if (flagout | flag2na) {
    if (!"flagout" %in% assayDataElementNames (object)) warning ("object does not have a flagout assayData element.")
  }
  
  ##taking out non "BiologicalFeature"
  if (nonbiologicalout | nonbiological2na) {
    if (!"BiologicalFeature" %in% colnames (fData(object))) {
      stop ("object does not have a BiologicalFeature column in pData. Nonbiological features cannot be excluded")
    }
    ##     if (verbose) {
    ##       message (paste (sum (!fData(object)[,"BiologicalFeature"]), "non biological features excluded"))
    ##       object <- object[fData(object)[,"BiologicalFeature"],]
    ##     }
  } 
  
  
  ##BG CORRECTING --------------------------------------------------------------

  if (verbose) cat ("\n")
  
  if (method == "agilent") {
      if (verbose) message ("using Agilent gProcessedSignal as background corrected intensities")
      bgcorrected.mat <- assayDataElement (object, "gProcessedSignal")

  } else {

    if (verbose) message (paste ("background correction method:",  method))
    
    ## foreground and background
    fg <- foreg (object)
    bg <- backg (object)
    if (is.null (bg)) {
      bg <- 0
    }
    
    ##several corrections
    if (method == "none") {
      bgcorrected.mat <- fg
    }
    if (method == "subtract") {
      bgcorrected.mat <- fg - bg
    }
    if (method == "half") {
        bgcorrected.mat <- fg - bg/2
    }

    ##RMA corrections
    if (method %in% c("normexp", "rma", "normexpnobg", "rmanobg")) {
        
        ##require (limma)
        
        if (method == "rma") {
        ##method <- "rma"
        dif <- fg - bg
      }
      if (method == "normexp") {
        method <- "saddle" #estilo limma; backgroundCorrect para matrices
        dif <- fg - bg
      }
      if (method == "rmanobg") {
        method <- "rma"
        dif <- fg #- bg
      }
      if (method == "normexpnobg") {
        method <- "saddle" #estilo limma; backgroundCorrect para matrices
        dif <- fg #- bg
      }
      
      bgcorrected.mat <- matrix (data = NA, nrow = nrow (dif), ncol = ncol (dif), dimnames = dimnames (dif))
      for (i in 1:ncol (dif)) {

        if (verbose) {
          cat ("\n")
          message (paste ("processing sample", i))
        }
        
        ##juego de flags
        if (flagout | nonbiologicalout) {
          
          excludedByFlag <- flagout & flags (object)[,i] #nos da la dimension necesaria
          
          if (nonbiologicalout) {#existe la columna "BiologicalFeature"
            excludedNonBiological <- nonbiologicalout & !fData(object)[,"BiologicalFeature"] #nonbiologicalout & ES REDUNDANTE
          } else {
            excludedNonBiological <- rep (FALSE, times = nrow (object))
          }

          touse <- !excludedByFlag & !excludedNonBiological
          
          if (verbose) {
            message (paste (sum (!touse), "spots not used to fit the normalization model"))
            print (table (excludedByFlag, excludedNonBiological))
          }
          
        } else {
          touse <- TRUE
        }

        ##esto es lo mas importante
        ne <- normexp.fit (dif[touse,i], method = method)
        bgcorrected.mat[,i] <- normexp.signal (ne$par, x = dif[,i])
        
      }
    }
    ##     if (method == "rma") {#rma eslilo limma
    ##       message (paste ("background correction method:",  method))
    ##       require (affy)
    ##       dif <- fg - bg
    ##       bgcorrected.mat <- apply (dif, 2, bg.adjust)
    ##     }
  }
  
  ##log ------------------------------------------------------------------------

  if (logit) {
    bgcorrected.mat <- log2 (bgcorrected.mat)
    if (verbose) message ("\nlog2 transformation done after background correction")
  }

  
  ##completamos NAs ------------------------------------------------------------
  
  if (flag2na) {
    if (verbose) message ("\nsetting flagged spots to NA")
    bgcorrected.mat[as.logical (flags (object))] <- NA
  }
  
  if (nonbiological2na) {
    if (verbose) message ("\nsetting Non Biological Features to NA")
    bgcorrected.mat[!fData(object)[,"BiologicalFeature"]] <- NA
  }

  
  ##keep: eliminate some assayDataElements -------------------------------------
  if (!keep) {
    elementos <- setdiff (assayDataElementNames (object), "exprs")
    if (verbose) message (paste ( c("\nElements: {", elementos, "} will be removed from the assayData"), collapse = " "))
    for (e in elementos) {
      object <- assayDataElementReplace (object, e, NULL)
    }
  }
  
  
  ##guardamos el "ExpressionSet" -----------------------------------------------

  object <- assayDataElementReplace (object, "exprs", bgcorrected.mat)
  
  ##chequeamos que el objeto de clase "ExpressionSet" es valido
  if (!validObject (object)) stop ("ExpressionSet cannot be modified")
  
  ##SALIDA
  return (object)
}


####################################################################################################################
####################################################################################################################


## library (agilent,  lib.loc = "/home/david/programas/mislibrerias/agilent/agilent_instalacion_local"); packageDescription ("agilent", fields = "Version") #
## setwd (file.path (system.file ("exampledata", package = "agilent")))
## dir ()

## ##datos <- readAgilent (background.column = NULL, other.columns = NULL, feature.columns = NULL, verbose = FALSE)
## ##object <- datos

## ##object <- readAgilent ()

## datos <- readGPR ()
## assayDataElementNames (datos)
## fData (datos)[,"BiologicalFeature"] <- fData (datos)[,"ControlType"] == "ignore"

## validObject (datos)

## datosA <- backgroundCorrectAgilent (datos)
## datosB <- backgroundCorrectAgilent (datos, verbose = FALSE)

## datosC <- backgroundCorrectAgilent (datos, method = "normexp")
## datosD <- backgroundCorrectAgilent (datos, method = "normexp", logit = FALSE)

## datosE <- backgroundCorrectAgilent (datos, method = "normexp", flagout = FALSE)
## datosF <- backgroundCorrectAgilent (datos, method = "normexp", flagout = FALSE, nonbiologicalout = FALSE)

## datosG <- backgroundCorrectAgilent (datos, method = "normexp", flagout = FALSE, nonbiologicalout = FALSE, flag2na = TRUE)
## datosH <- backgroundCorrectAgilent (datos, method = "normexp", flagout = FALSE, nonbiologicalout = FALSE, flag2na = TRUE, nonbiological2na = TRUE)
## datosI <- backgroundCorrectAgilent (datos, method = "normexp", flagout = FALSE, nonbiologicalout = FALSE, nonbiological2na = TRUE)
## datosJ <- backgroundCorrectAgilent (datos, method = "normexp", flagout = TRUE,  nonbiologicalout = TRUE, flag2na = TRUE, nonbiological2na = TRUE)

## methods <- c ("agilent", "none", "subtract", "half", "normexp", "rma")
## datosM <- backgroundCorrectAgilent (datos, method = "agilent")
## datosN <- backgroundCorrectAgilent (datos, method = "none")
## datosO <- backgroundCorrectAgilent (datos, method = "subtract")
## datosP <- backgroundCorrectAgilent (datos, method = "half")
## datosQ <- backgroundCorrectAgilent (datos, method = "rma")

## datosZ <- backgroundCorrectAgilent (datos, keep = FALSE)

## i <- 1

## assayDataElementNames (datosB)
## plot (exprs (datosF)[,i], exprs (datosC)[,i])

## table (exprs (datosB) == exprs (datosC))

## plot (exprs (datosB)[,i],  exprs (datosC)[,i])
## plot (exprs (datosD)[,i],  exprs (datosC)[,i])

## summary (exprs (datosB))
## summary (exprs (datosC))
## summary (exprs (datosD))
## summary (exprs (datosE))
## summary (exprs (datosF))
## summary (exprs (datosG))
## summary (exprs (datosH))
## summary (exprs (datosJ))

## summary (exprs (datosM))
## summary (exprs (datosN))
## summary (exprs (datosO))
## summary (exprs (datosP))
## summary (exprs (datosQ))

## summary (exprs (datosZ))

## apply (exprs (datosJ), 2,function (x) sum (is.na (x)))

## table (is.na ((exprs (datosG))) == flags (datos))

## exprs (datosJ)[1:3,]

## assayDataElementNames (datosA)
## assayDataElementNames (datosZ)
