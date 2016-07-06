##deriveAgilentFeatures.r
##2009-08-23 dmontaner@cipf.es
##2010-10-29 dmontaner@cipf.es
##2013-11-29 dmontaner@cipf.es

## FALTA PONER CHEQUEO DE QUE ESTAN TODAS LAS COLUMNAS QUE SE NECESITAN
## c ("SubTypeName", "SubTypeMask", "ControlType", "SpikeInlogRC", "SystematicName")

##' @name deriveAgilentFeatures
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @keywords agilent features
##' @seealso \code{\link{readAgilent}} \code{\link{readGPR}}
##'
##' @title Derive Agilent Feature Data.
##'
##' @description Describes Agilent Feature data.
##' 
##' @details Derives some feature data for the spots of the array.
##' Agilent features "SubTypeName" and "ControlType" are combined into a new "FeatureType".
##' Notice that Agilent features "SubTypeName" and "SubTypeMask" are equivalent.
##' Only spots tagged as "BiologicalFeature" in the "FeatureType" are suitable for analyzing.
##' Spots tagged as "E1A" in the "FeatureType" are SpikeIns.
##' log.relative.conc is a named vector describing
##' 1-color Agilent SpikeIns Signal Statistics Log (Relative Conc.)
##' Concentrations form Agilent Design ID 014850 are assumed by default.
##' 
##' @param data a data.frame of feature Data obtained from readAgilent
##' @param log.relative.conc named vector of SpikeIns log relative concentration. 
##' @param verbose verbose
##' 
##' @return A data.frame with derived features.
##'
##' @export

deriveAgilentFeatures <- function (data, log.relative.conc = NULL, verbose = TRUE) {

  ##checking that the necessary columns are in the data
  neededcol <- c ("SubTypeName", "SubTypeMask", "ControlType", "SystematicName")
  if (!all (neededcol %in% colnames (data))) {
    stop (paste (c("columns", setdiff (neededcol, colnames (data)), "are not found in data"), collapse = " "))
  }
          
  ##----------------------------------------------------------------------------
  ##SpikeIns - ESTO DEBERIACAMBIAR CON LOS TIPOS DE CHIPS
  ##1-color Agilent SpikeIns Signal Statistics
  ##Log (Relative Conc.)

  if (is.null (log.relative.conc)) {
    log.relative.conc <- c (E1A_r60_3    = 0.3,
                            E1A_r60_a104 = 1.3,
                            E1A_r60_a107 = 2.3,
                            E1A_r60_a135 = 3.3,
                            E1A_r60_a20  = 3.83,
                            E1A_r60_a22  = 4.3,
                            E1A_r60_a97  = 4.82,
                            E1A_r60_n11  = 5.3,
                            E1A_r60_n9   = 5.82,
                            E1A_r60_1    = 6.3)
    if (verbose) message ("Understanding SpikeIns as described for Agilent Design ID 014850")
  }
  ##----------------------------------------------------------------------------
  ##"SpikeIns as described for Agilent Design ID 014850"
  
  ##creating "FeatureType" describing each spot
  ##summarizes "SubTypeName" and "ControlType" variables from the original data
  ##notice that "SubTypeName" is equivalent to "SubTypeMask"
  feature.type <- data[,"SubTypeName"]
  ##feature.type[is.na (feature.type)] <- "BiologicalFeature"
  feature.type[is.na (feature.type) & data[,"ControlType"] ==  0] <- "BiologicalFeature"
  feature.type[data[,"SubTypeMask"] == 66 & data[,"ControlType"] ==  1] <- "StructuralPositive"
  feature.type[data[,"SubTypeMask"] == 66 & data[,"ControlType"] == -1] <- "StructuralNegative"
  ##feature.type[data[,"SubTypeMask"] == 8196 & data[,"ControlType"] ==  1] <- "Spikein"
  feature.type[is.na (feature.type)] <- "unknown"
  
  ##creating indicator variables for each "FeatureType"
  boolean <- sapply (unique (feature.type), function (u) feature.type == u, USE.NAMES = TRUE)
  
  ##combining into a data.frame
  derivedvars <- cbind (data.frame (FeatureType = feature.type, stringsAsFactors = FALSE),
                        as.data.frame (boolean, stringsAsFactors = FALSE))
  
  ##adding SpikeIns concentration according to "SystematicName"
  derivedvars[,"SpikeInlogRC"] <- NA
  derivedvars[feature.type == "E1A", "SpikeInlogRC"] <-
    log.relative.conc[data[feature.type == "E1A", "SystematicName"]]

  if (verbose) {
    cat ("\nSpikeIn Log Relative Conc.\n")
    print (table (derivedvars[,"SpikeInlogRC"], exclude = NULL))
  }
  
  ##identificamos los "ProbeName" que estan duplicados
  derivedvars[,"ProbeNameTimes"]  <- table (data[,"ProbeName"])[data[,"ProbeName"]]
  derivedvars[,"uniqueProbeName"] <- derivedvars[,"ProbeNameTimes"] == 1
  derivedvars[,"hasDescription"]  <- !is.na (data[,"Description"]) & data[,"Description"] != "Unknown"
 
  ##SALIDA
  return (derivedvars)
}

## setwd ("/home/david/datos/2009/amparo_galan/data_processed")
## load ("raw_data.RData")
## ls ()

## out <- deriveAgilentFeatures (fData (datos.raw))

## table (out == fData (datos.raw)[,colnames (out)])

## f <- function (x) {
##   print (is.null(x))
## }
## f()
## f(3)
## f(NULL)
