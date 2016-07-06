##readAgilent.r
##2009-10-17 dmontaner@cipf.es
##2010-05-18 dmontaner@cipf.es
##2013-11-29 dmontaner@cipf.es

##FALTA ACABAR LA AYUDA
##INCLUIR UNA VARIABLE PATH

##' @name readAgilent
##' @author David Montaner \email{dmontaner@@cipf.es}
##'
##' @keywords read agilent
##' @seealso \code{\link{readAgilentHeader}} \code{\link{readGPR}}
##'
##' @title Read Agilent.
##' 
##' @description Reads Agilent raw data txt exported files.
##' 
##' @details exprs in assayData is set to NA
##'
##' sampleinfo is used rather than files or samplenames
##'
##' sampleinfo must have columns named 'filename' and 'samplename'
##'
##' If background.columncolumn = NULL no B matrix is set in assayData.
##' If all feature.columns are the same a single column is passed to the featureData.
##' If not all of them are equal, the information stays in the assayData
##' 
##' Default feature.columns are: "ProbeName", "GeneName", "SystematicName", "FeatureNum",
##' "Row", "Col", "SubTypeMask", "SubTypeName", "ProbeUID", "ControlType", "Description"
##'
##' Extra columns read by default are: "gProcessedSignal" and "IsManualFlag".
##' See the \code{other.columns} parameter.
##' 
##'
##' @param files names of the files to be read.
##' @param samplenames names to be given to the arrays.
##' @param sampleinfo data frame with columns 'filename', 'samplename' and any other sample information.
##' @param foreground.column column name of the foreground intensities.
##' @param background.column column name of the background intensities.
##' @param other.columns some other columns to be read.
##' @param feature.columns some other columns to be read as featureData.
##' @param derivefeatures should additional features be computed.
##' @param verbose verbose
##' @param \dots further arguments to be passed to deriveAgilentFeatures.
##' 
##' @return An ExpressionSet containing foreground (F) and background (B) matrices within the assayData
##' 
##' @examples
##' setwd (file.path (system.file ("exampledata", package = "agilent")))
##' dir ()
##' 
##' ra <- readAgilent ()
##' pData (ra)
##' fData (ra)[1:5,]
##' assayDataElementNames (ra)
##' 
##' sinfo <- cbind (filename = c("agilent2.txt", "agilent1.txt"),
##'                 samplename = c("first", "second"), more = 1:2)
##' ra <- readAgilent (sampleinfo = sinfo)
##' pData (ra)
##' 
##' ra <- readAgilent (background.column = NULL, other.columns = NULL, feature.columns = NULL, verbose = FALSE)
##' pData (ra)
##' fData (ra)[1:5,]
##' assayDataElementNames (ra)
##'
## @import Biobase
##' 
##' @export

## files = dir (pattern = ".txt")
## samplenames = NULL
## sampleinfo = NULL
## foreground.column = "gMeanSignal"
## background.column = "gBGMedianSignal"
## other.columns = NULL
## feature.columns = NULL
## derivefeatures = TRUE
## verbose = TRUE

readAgilent <- function (files = character (0), samplenames = NULL, sampleinfo = NULL,
                         foreground.column = "gMeanSignal", background.column = "gBGMedianSignal",
                         other.columns = c ("gProcessedSignal", "IsManualFlag"),
                         feature.columns = character(0),
                         derivefeatures = TRUE, verbose = TRUE, ...) {
  
  ##require (Biobase)
  if (verbose) cat ("\n") #nice output


  ##default files --------------------------------------------------------------
  if (length (files) == 0) {
    files <- dir (pattern = ".txt", ignore.case = TRUE)
  }
  
  ##default feature.columns ----------------------------------------------------
  ##meterlo en la definicion de la funcion daba problemas al roxygen
  if (length (feature.columns) == 0) {
    feature.columns <- c ("ProbeName", "GeneName", "SystematicName", "FeatureNum",
                          "Row", "Col", "SubTypeMask", "SubTypeName",
                          "chr_coord", "probe_mappings",
                          "ProbeUID", "ControlType", "Description")
  }

  
  ##organizing filenames and samplenames ---------------------------------------

  if (is.null (sampleinfo)) {
    if (is.null (samplenames)) {
      samplenames <- sub (pattern = ".txt", replacement = "", files, ignore.case = TRUE)
    }
    sampleinfo <- as.data.frame (cbind (filename = files, samplename = samplenames),
                                 stringsAsFactors = FALSE)
  } else {
    if (verbose) {
      message ("\nsampleinfo is used rather than files or samplenames\n")
    }
    if (!all (c ("filename", "samplename") %in% colnames (sampleinfo))){
      stop ("sampleinfo must have columns filename and samplename")
    }
    files       <- sampleinfo[,"filename"]
    samplenames <- sampleinfo[,"samplename"]
  }
  rownames (sampleinfo) <- sampleinfo[,"samplename"]
  sampleinfo <- as.data.frame (sampleinfo, stringsAsFactors = FALSE)

  
  ##reading headers --------------------------------------------------------------
  headers <- readAgilentHeaderS (files = files, samplenames = samplenames,
                                 includeFileNames = FALSE, verbose = verbose)
  if (verbose) cat ("\n") #nice output
  
  ##checking that columns are of the same type across microarrays
  if (is.data.frame (headers$fieldType)) stop ("files seem to have different column layout")

  
  ##checking some parameters ---------------------------------------------------
  
  if (length (foreground.column) > 1) {
    warning ("foreground.column has more than one element; just the first one will be used")
    foreground.column <- foreground.column[1]
  }
  if (!foreground.column %in% names (headers$fieldType)) {
    stop (paste (foreground.column, "does not seem to be a column in your data"))
  }

  if (!is.null (background.column)) {
    if (length (background.column) > 1) {
      warning ("background.column has more than one element; just the first one will be used")
      background.column <- background.column[1]
    }
    if (!background.column %in% names (headers$fieldType)) {
      stop (paste (background.column, "does not seem to be a column in your data"))
    }
  }
  
  bad.other.columns <- setdiff (other.columns, names (headers$fieldType))
  if (length (bad.other.columns) > 0) {
    warning (paste (c ("columns", bad.other.columns, "are not found in your data"), collapse = " "))
  }
  other.columns <- intersect (names (headers$fieldType), other.columns)
  ##OBS: THE intersect HAS TO BE STRICTLY IN THIS ORDER
  ##     IF WE WANT A null VALUE TO YIELD ALSO A null INTERSECTION
  
  bad.feature.columns <- setdiff (feature.columns, names (headers$fieldType))
  if (length (bad.feature.columns) > 0) {
    warning (paste (c ("columns", bad.feature.columns, "are not found in your data"), collapse = " "))
  }
  feature.columns <- intersect (names (headers$fieldType), feature.columns) 
  ##OBS: THE intersect HAS TO BE STRICTLY IN THIS ORDER
  ##     IF WE WANT A null VALUE TO YIELD ALSO A null INTERSECTION

  
  ## READING -------------------------------------------------------------------
  
  ##defining columns to be read
  cols2read <- unique (c(foreground.column, background.column, other.columns, feature.columns))
  
  ##creating assayData for the "ExpressionSet"
  assaydata <- assayDataNew (storage.mode = "environment") #"lockedEnvironment" no se puede sobreescribir
  for (v in cols2read) assign (v, value = NULL, envir = assaydata)

  ##colum class. Setting as NULL columns that should not be read
  col.class <- c (primerafila = "NULL", headers$fieldType)
  col.class[!names (col.class) %in% cols2read] <- "NULL"
  col.class <- sub ("logical", "integer", col.class)

  
  ##reading files --------------------------------------------------------------

  if (verbose) {
    message ("Reading:")
    print (cbind (files, samplenames))
  }
  
  for (fichero in files) {

    if (verbose) {
      cat ("\n")
      message (paste ("Reading:", fichero))
    }
    
    read.features <- read.table (file = fichero, header = TRUE, sep = "\t", quote = "",
                                 na.strings = "", colClasses = col.class,
                                 skip = 9, comment.char = "")

    ##defining logical variables
    was.logical <- headers$fieldType[colnames (read.features)] == "logical"
    for (col in colnames (read.features)[was.logical]) { 
      read.features[,col] <- as.logical (read.features[,col])
    }

    if (verbose) message (paste ("dim", nrow (read.features), ncol (read.features)))
  
    ##merging data into the assaydata environment
    for (v in cols2read) assaydata[[v]] <- cbind (assaydata[[v]], read.features[,v])
  }

  
  ##parsing feature data -------------------------------------------------------

  if (is.null (feature.columns)) {
    featuredata <- as.data.frame (matrix (nrow = nrow (read.features), ncol = 0))
  } else {
    if (verbose)  message ("\nparsing feature.columns")
    featuredata <- list ()
    for (f in feature.columns) {
      if (verbose) print (f)
      ml <- max (apply (assaydata[[f]], 1, function (x) length (unique (x))))
      if (ml == 1) {
        mylist <- list (assaydata[[f]][,1])
        names (mylist) <- f
        featuredata <- c(featuredata, mylist)
        rm (list = f, envir = assaydata)  ##remove matrix from the assaydata environment
      } else {
        warning (paste (f, "has not the same values across all microarrays. Values stay in the assayData"))
      }
    }
    featuredata <- as.data.frame (featuredata, stringsAsFactors = FALSE)
  }


  ## derive features some other features ----------------------------------------

  ##checking that all columns needed to derive features are available
  neededcol <- c ("SubTypeName", "SubTypeMask", "ControlType", "SystematicName")
  
  if (derivefeatures) {
    if (!all (neededcol %in% colnames (featuredata))) {
      warning ("further features could not be derived")
      if (verbose) {
        cat ("\n") #nice output
        message ("further features could not be derived")
        cat ("\n") #nice output
        message (paste (c("columns", setdiff (neededcol, colnames (featuredata)),
                          "are needed to derive further features"), collapse = " "))
        cat ("\n") #nice output
      }
      derivefeatures <- FALSE
    }
  }

  if (derivefeatures) {
    if (verbose) {
      cat ("\n") #output formatting
      message ("deriving Agilent features")
    }
    featuredata <- cbind (featuredata, deriveAgilentFeatures (featuredata, verbose = verbose, ...))
  } else {
    featuredata[,"BiologicalFeature"] <- featuredata[,"ControlType"] == 0
  }
    
  
  ##setting Foreground and Background ------------------------------------------
  
  assaydata[["F"]] <- assaydata[[foreground.column]]
  rm (list = foreground.column, envir = assaydata)

  if (!is.null (background.column)) {
    assaydata[["B"]] <- assaydata[[background.column]]
    rm (list = background.column, envir = assaydata)
  }
  
  ##building exprs matrix ------------------------------------------------------
  assaydata[["exprs"]] <-  matrix (data = NA, nrow = nrow (read.features), ncol = length (samplenames))
  ##   if (create.exprs) {
  ##     if (verbose) {
  ##       cat ("\n") #output formatting
  ##       message ("exprs matrix is created within the assayData")
  ##     }
  ##     assaydata[["exprs"]] <-  matrix (data = NA, nrow = nrow (read.features), ncol = length (samplenames))
  ##   }
  
  
  ##setting col and row names ----------------------------------------------------------

  for (v in ls (assaydata)) {
    colnames (assaydata[[v]]) <- samplenames
    rownames (assaydata[[v]]) <- rownames (featuredata)
  }

  
  ##closing the environment ----------------------------------------------------

  storageMode (assaydata) <- "lockedEnvironment"

  
  ##building "ExpressionSet" -----------------------------------------------------

  expressionset <- new ("ExpressionSet", assayData = assaydata,
      phenoData = new ("AnnotatedDataFrame", data = sampleinfo,
        varMetadata = data.frame (labelDescription = rep ("from sampleinfo", times = ncol (sampleinfo)),
          row.names = colnames (sampleinfo), stringsAsFactors = FALSE)),
      featureData = new ("AnnotatedDataFrame", data = featuredata,
        varMetadata = data.frame (labelDescription = rep ("read from file", times = ncol (featuredata)),
          row.names = colnames (featuredata), stringsAsFactors = FALSE)))
  
  ##checking the "ExpressionSet"
  if (!validObject (expressionset)) warning ("not valid ExpressionSet created")


 ##OUT
  return (expressionset)
}

################################################################################
################################################################################

## library (agilent, lib.loc = "/home/david/programas/mislibrerias/agilent/agilent_instalacion_local"); packageDescription ("agilent", fields = "Version") #
## searchpaths ()
## help (package = agilent)
## setwd (file.path (system.file ("exampledata", package = "agilent")))
## dir ()

## ra <- readAgilent ()
## pData (ra)
## fData (ra)[1:5,]
## assayDataElementNames (ra)

## ra <- readAgilent (derivefeatures = FALSE) ##OK
## ra <- readAgilent (feature.columns = c("ProbeName", "ProbeUID")) ##OK

## ra <- readAgilent (files = c("agilent2.txt", "agilent1.txt"), samplenames = c("dos", "uno")) #OK
## pData (ra)
## sapply (pData (ra), class)


## sinfo <- cbind (filename = c("agilent2.txt", "agilent1.txt"),
##                 samplename = c("first", "second"), more = 1:2)
## ra <- readAgilent (sampleinfo = sinfo)
## pData (ra)
## sapply (pData (ra), class)

## ra <- readAgilent (feature.columns = c("ProbeName", "ProbeUID"), derivefeatures = FALSE) ##OK pero

## ra <- readAgilent (foreground.column = "gProcessedSignal", feature.columns = c("ProbeName", "ProbeUID", "gMeanSignal"), derivefeatures = FALSE) ##OK pero
## fData (ra)[1:5,]
## assayDataElementNames (ra)

## ra <- readAgilent (background.column = NULL, other.columns = NULL, feature.columns = NULL)
## pData (ra)
## fData (ra)[1:5,]
## assayDataElementNames (ra)

## readAgilent (verbose = FALSE)

## debug (readAgilent)
## undebug (readAgilent)
