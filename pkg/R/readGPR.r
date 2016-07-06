##readGPR.r
##2009-10-19 dmontaner@cipf.es
##2010-10-24 dmontaner@cipf.es
##2013-11-29 dmontaner@cipf.es

##reading GenePix raw data files
##FALTA ACABAR LA AYUDA
##INCLUIR UNA VARIABLE PATH



##' @name readGPR
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @keywords Read GPR.
##' @seealso \code{\link{readAgilent}} \code{\link{readGPRHeaderS}}
##' 
##' @title Average duplicated rows
##'
##' @description Reads GPR raw data txt exported files.
##' 
##' @details exprs in assayData is set to NA.
##'
##' sampleinfo is used rather than files or samplenames.
##'
##' sampleinfo must have columns filename and samplename.
##'
##' If background.column = NULL no B matrix is set in assayData
##' If derivefeatures is TRUE "BiologicalFeature" is added to he fData indicating when ControlType is "FALSE".
##'
##' If all feature.columns are the same a single column is passed to the featureData.
##' If not all of them are equal, the information stays in the assayData
##' 
##' default feature.columns are "Block", "Column", "Row", "Name", "ID" and "ControlType"
##' For consistency with the original Agilent text format, if "Column" is among
##' feature.columns a final column named just "Col" is in the pData object.
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
##' @param \dots further arguments NOT IN USE
##' 
##' @return An ExpressionSet containing foreground (F) and background (B) matrices within the assayData.
##' 
##' @examples
##' setwd (file.path (system.file ("exampledata", package = "agilent")))
##' dir ()
##' 
##' rg <- readGPR ()
##' pData (rg)
##' fData (rg)[1:5,]
##' assayDataElementNames (rg)
##' table (assayData (rg)[["Flags"]][,1] == assayData (rg)[["Flags"]][,2])
##' 
##' sinfo <- cbind (filename = c("GPR2.gpr", "GPR1.gpr"),
##'                   samplename = c("second", "first"), more = 1:2)
##' rg <- readGPR (sampleinfo = sinfo)
##' pData (rg)
##' 
##' rg <- readGPR (background.column = NULL, other.columns = NULL, feature.columns = NULL, verbose = FALSE)
##' pData (rg)
##' fData (rg)
##' assayDataElementNames (rg)
##' exprs (rg) [1:3,]
##' foreg (rg) [1:3,]
##'
## @import Biobase
##' 
##' @export

################################################################################
################################################################################

## rm (list = ls ())
## library (agilent,  lib.loc = "/home/david/programas/mislibrerias/agilent/agilent_instalacion_local"); packageDescription ("agilent", fields = "Version") #
## searchpaths ()
## setwd (file.path (system.file ("exampledata", package = "agilent")))

## files = character (0)
## samplenames = NULL
## sampleinfo = NULL
## foreground.column = "F532 Mean"
## background.column = "B532 Median"
## background.column = NULL
## other.columns = "Flags"
## other.columns = NULL
## feature.columns = c("Block", "Column", "Row", "Name", "ID", "ControlType") ##, "X", "Y"
## feature.columns = NULL
## derivefeatures = TRUE
## verbose = TRUE

readGPR <- function (files = character (0), samplenames = NULL, sampleinfo = NULL,
                     foreground.column = "F532 Mean", background.column = "B532 Median",
                     other.columns = c ("Flags"),
                     feature.columns = c("Block", "Column", "Row", "Name", "ID", "ControlType"),
                     derivefeatures = TRUE, verbose = TRUE, ...) {#2009-10-19
  
  ##require (Biobase)
  if (verbose) cat ("\n") #nice output
  
  
  ##default files --------------------------------------------------------------
  if (length (files) == 0) {
    files <- dir (pattern = ".gpr", ignore.case = TRUE)
  }
  
  ##default feature.columns ----------------------------------------------------
  ##meterlo en la definicion de la funcion daba problemas al roxygen
##   if (length (feature.columns) == 0) {
##     feature.columns <- c ("ProbeName", "GeneName", "SystematicName", "FeatureNum",
##                           "Row", "Col", "SubTypeMask", "SubTypeName", #"chr_coord",
##                           "ProbeUID", "ControlType", "Description")
##   }

  
  ##organizing filenames and samplenames ---------------------------------------

  if (is.null (sampleinfo)) {
    if (is.null (samplenames)) {
      samplenames <- sub (pattern = ".gpr", replacement = "", files, ignore.case = TRUE)
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

  
  ##reading headers ------------------------------------------------------------
  headers <- readGPRHeaderS (files = files, samplenames = samplenames,
                             includeFileNames = TRUE, verbose = verbose)
  if (verbose) cat ("\n") #nice output
  
  ##checking that columns are of the same type across microarrays
  ##if (length (unique (headers$params$NHeaderRecords)) > 1) stop ("files seem to have different column layout")
  if (length (unique (headers$params$NHeaderRecords)) > 1) warning ("files seem to have different header layout")
  
  
  ##checking some parameters ---------------------------------------------------
  
  if (length (foreground.column) > 1) {
    warning ("foreground.column has more than one element; just the first one will be used")
    foreground.column <- foreground.column[1]
  }
  if (!foreground.column %in% headers$columns) {
    stop (paste (foreground.column, "does not seem to be a column in your data"))
  }

  if (!is.null (background.column)) {
    if (length (background.column) > 1) {
      warning ("background.column has more than one element; just the first one will be used")
      background.column <- background.column[1]
    }
    if (!background.column %in% headers$columns) {
      stop (paste (background.column, "does not seem to be a column in your data"))
    }
  }
  
  bad.other.columns <- setdiff (other.columns, headers$columns)
  if (length (bad.other.columns) > 0) {
    warning (paste (c ("columns", bad.other.columns, "are not found in your data"), collapse = " "))
  }
  other.columns <- intersect (headers$columns, other.columns)
  ##OBS: THE intersect HAS TO BE STRICTLY IN THIS ORDER
  ##     IF WE WANT A null VALUE TO YIELD ALSO A null INTERSECTION
  
  bad.feature.columns <- setdiff (feature.columns, headers$columns)
  if (length (bad.feature.columns) > 0) {
    warning (paste (c ("columns", bad.feature.columns, "are not found in your data"), collapse = " "))
  }
  feature.columns <- intersect (headers$columns, feature.columns) 
  ##OBS: THE intersect HAS TO BE STRICTLY IN THIS ORDER
  ##     IF WE WANT A null VALUE TO YIELD ALSO A null INTERSECTION
    
  
  ## READING -------------------------------------------------------------------
  
  ##defining columns to be read
  foreground.column <- make.names (foreground.column, unique = TRUE)
  background.column <- make.names (background.column, unique = TRUE)
  other.columns     <- make.names (other.columns, unique = TRUE)
  feature.columns   <- make.names (feature.columns, unique = TRUE)
  
  cols2read <- unique (c(foreground.column, background.column, other.columns, feature.columns))

  
  ##creating assayData for the "ExpressionSet"
  assaydata <- assayDataNew (storage.mode = "environment") #"lockedEnvironment" no se puede sobreescribir
  for (v in cols2read) assign (v, value = NULL, envir = assaydata)

##   ##colum class. Setting as NULL columns that should not be read
##   col.class <- c (primerafila = "NULL", headers$fieldType)
##   col.class[!names (col.class) %in% cols2read] <- "NULL"
##   col.class <- sub ("logical", "integer", col.class)

  
  ##reading files --------------------------------------------------------------
  ##NOTE read.table reads the FALSE values in the ControlType column as LowerCase character "false"
  
  if (verbose) {
    message ("Reading:")
    print (cbind (files, samplenames))
  }
  
  for (fichero in files) {

    if (verbose) {
      cat ("\n")
      message (paste ("Reading:", fichero))
    }
    
    myskip <- headers$params[headers$params$fileName == fichero, "NHeaderRecords"]
    
    read.features <- read.table (file = fichero, header = TRUE, sep = "\t",
                                 quote = "\"", as.is = TRUE, skip = myskip)
        
    ##merging data into the assaydata environment
    for (v in cols2read) assaydata[[v]] <- cbind (assaydata[[v]], read.features[,v])

    ##printing the dimension of the read data
    if (verbose)  {
      read.features <- as.data.frame (read.features[,cols2read]) #problema cuando solo hay una variable se convierte en un vector
      message (paste ("dim", nrow (read.features), ncol (read.features)))
    }
  }

  ##flags ----------------------------------------------------------------------
  ##esto habria que repasarlo
  if (derivefeatures) {
    if ("Flags" %in% assayDataElementNames (assaydata)) {
      assaydata[["flagout"]] <- assaydata[["Flags"]] < 0
    }
  }
  
  ##parsing feature data -------------------------------------------------------
  
  if (length (feature.columns) == 0) {
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
  
  ##change the name "Column" to "Col" as in the Agilent scanner layout
  colnames (featuredata)[colnames (featuredata) == "Column"] <- "Col"

  
  ## derive features some more features ----------------------------------------
  if (derivefeatures) {
    if (verbose) {
      cat ("\n") #output formatting
      message ("deriving GPR features")
    }
    if ("ControlType" %in% colnames (featuredata)) {
      featuredata <- cbind (featuredata, BiologicalFeature = tolower (featuredata$ControlType) == "false")
      if (verbose) {
        message ("BiologicalFeature filed created")
      }
    } else {
      if (verbose) {
        message ("ControlType column was not found. BiologicalFeature filed could not be created")
      }
    }
  }
  
##   ##checking that all columns needed to derive features are available
##   neededcol <- c ("SubTypeName", "SubTypeMask", "ControlType", "SystematicName")
  
##   if (derivefeatures) {
##     if (!all (neededcol %in% colnames (featuredata))) {
##       warning ("further features could not be derived")
##       if (verbose) {
##         cat ("\n") #nice output
##         message ("further features could not be derived")
##         cat ("\n") #nice output
##         message (paste (c("columns", setdiff (neededcol, colnames (data)),
##                           "are needed to derive further features"), collapse = " "))
##         cat ("\n") #nice output
##       }
##       derivefeatures <- FALSE
##     }
##   }

##   if (derivefeatures) {
##     if (verbose) {
##       cat ("\n") #output formatting
##       message ("deriving GPR features")
##     }
##     featuredata <- cbind (featuredata, deriveGPRFeatures (featuredata, ...))
##   }
  
  
  ##setting Foreground and Background ------------------------------------------
  
  assaydata[["F"]] <- assaydata[[foreground.column]]
  rm (list = foreground.column, envir = assaydata)
  
  if (length (background.column) == 1) {
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

## library (agilent,  lib.loc = "/home/david/programas/mislibrerias/agilent/agilent_instalacion_local"); packageDescription ("agilent", fields = "Version") #

## debug (readGPR)
## undebug (readGPR)
