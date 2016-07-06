##readAgilentHeaderS.r
##2009-08-22 dmontaner@cipf.es
##2009-11-21 dmontaner@cipf.es
##2013-11-29 dmontaner@cipf.es


##' @name readAgilentHeaderS
##' @author David Montaner \email{dmontaner@@cipf.es}
##'
##' @keywords read agilent header
##' @seealso \code{\link{readAgilentHeader}} \code{\link{readAgilent}} \code{\link{readGPR}}
##'
##' @title Read Agilent headers.
##' 
##' @description Reads Agilent header from one or several txt exported files.
##' 
##' @details FEPARAMS and STATS from header of all arrays are stored in a data.frames.
##' FEATURE TYPE is stored in a character vector (named fieldType) when it is the same for all samples.
##' If it is not the same across arrays (or if includeFileNames = TRUE) a data.frame is used an a warning message issued.
##' Agilent \strong{field types} are coded in \R standards.
##' 
##' @param files names of the files which header is to be read.
##' @param samplenames names to be given to the arrays.
##' @param verbose verbose
##' @param includeFileNames should the filename be included into datasets?
##' 
##' @return A list of 3 components containing the FEPARAMS, STATS, and FEATURE TYPE FROM the header of the Agilent files.
##' 
##' @examples
##' setwd (file.path (system.file ("exampledata", package = "agilent")))
##' dir ()
##'
##' ra <- readAgilentHeaderS ()
##' names (ra)
##' dim (ra$feparams)
##' ra$feparams[,1:3]
##' class (ra[[3]])
##' ra[[3]]
##'
##' ra <- readAgilentHeaderS (files = c("agilent1.txt", "agilent2.txt"), samplenames = c ("first", "second"), includeFileNames = FALSE)
##' names (ra)
##' dim (ra$feparams)
##' ra$feparams[,1:3]
##' class (ra[[3]])
##'
##' @importFrom utils read.table
##' @export

readAgilentHeaderS <- function (files = dir (pattern = ".txt"), samplenames = sub (".txt", "", files),
                                verbose = TRUE, includeFileNames = FALSE) {#2009-11-21 dmontaner@cipf.es
    
  feparams  <- data.frame (INIT = NA)
  stats     <- data.frame (INIT = NA)
  fieldType <- data.frame (INIT = NA)
  
  for (fichero in files) {
    if (verbose) cat (paste ("reading header:", fichero), fill = TRUE)

    fileHeader <- readAgilentHeader (fichero)
    
    feparams  <- merge (feparams,  data.frame (fileHeader$feparams,  stringsAsFactors = FALSE), all = TRUE) # all = TRUE is important
    stats     <- merge (stats,     data.frame (fileHeader$stats,     stringsAsFactors = FALSE), all = TRUE)
    fieldType <- merge (fieldType, data.frame (fileHeader$fieldType, stringsAsFactors = FALSE), all = TRUE)
  }
      
  feparams  <- feparams [colnames (feparams)  != "INIT"]
  stats     <- stats    [colnames (stats)     != "INIT"]
  fieldType <- fieldType[colnames (fieldType) != "INIT"]

  if (!is.null (samplenames)) {
    rownames (feparams)  <- samplenames
    rownames (stats)     <- samplenames
    rownames (fieldType) <- samplenames
  }
  
  if (!includeFileNames) {
    feparams  <- feparams [colnames (feparams)  != "filename"]
    stats     <- stats    [colnames (stats)     != "filename"]
    fieldType <- fieldType[colnames (fieldType) != "filename"]
  }
  
  ##if fieldType is the same across arrays we return just a named vector
  if (all (apply (fieldType, 2, function (x) length (unique (x)) == 1))) {
    fieldType <- unlist (fieldType[1,])
  } else {
    warning ("Agilent FEATURES don not have the same TYPE across arrays")
    if (includeFileNames) {
      fieldType <- cbind (fileName = files, fieldType, stringsAsFactors = FALSE)
    }
  }
  
  ###RETURN
  return (list (feparams = feparams, stats = stats, fieldType = fieldType))
}

################################################################################

##' @name readAgilentHeader
##' @author David Montaner \email{dmontaner@@cipf.es}
##'
##' @keywords read agilent header
##' @seealso \code{\link{readAgilentHeaderS}} \code{\link{readAgilent}} \code{\link{readGPR}}
##'
##' @title Read Agilent header form a single file.
##' 
##' @description Reads Agilent header from a unique txt exported file.
##'
##' @details readAgilentHeaderS should be generally used.
##' Reads Agilent header from one or several txt exported files.
##' 
##' @param file name of the file which header is to be read.
##' 
##' @return A list of 3 components containing the FEPARAMS, STATS, and FEATURE TYPE FROM the header of the Agilent file.
##' Each element of the list is a list itself.
##' 
##' @examples
##' myfile <- file.path (system.file ("exampledata", package = "agilent"), "agilent1.txt")
##' myfile
##' readAgilentHeader (file = myfile)
##'
##' @export

readAgilentHeader <- function (file) {
  
  ##FEPARAMS
  feparams <- read.table (file = file, header = FALSE, sep = "\t", quote = "", colClasses = "character", nrows = 3)
  feparams <- feparams[,-1]
  ##
  feparams.mode  <- unlist (feparams[1,])
  feparams.names <- unlist (feparams[2,])
  ##
  feparams <- feparams[3,]
  colnames (feparams) <- feparams.names
  rownames (feparams) <- NULL
  ##
  for (i in which (feparams.mode == "text"))    mode (feparams[,i]) <- "character"  
  for (i in which (feparams.mode == "integer")) mode (feparams[,i]) <- "integer"
  for (i in which (feparams.mode == "float"))   mode (feparams[,i]) <- "numeric"
  for (i in which (feparams.mode == "boolean")) {mode (feparams[,i]) <- "integer"; mode (feparams[,i]) <- "logical"}
  ##
  feparams <- as.list (feparams)
  feparams["filename"] <- file

  ##STATS
  stats <- read.table (file = file, header = FALSE, sep = "\t", quote = "", colClasses = "character", nrows = 3, skip = 4)
  stats <- stats[,-1]
  ##
  stats.mode  <- unlist (stats[1,])
  stats.names <- unlist (stats[2,])
  ##
  stats <- stats[3,]
  colnames (stats) <- stats.names
  rownames (stats) <- NULL
  ##
  for (i in which (stats.mode == "text"))    mode (stats[,i]) <- "character"  
  for (i in which (stats.mode == "integer")) mode (stats[,i]) <- "integer"
  for (i in which (stats.mode == "float"))   mode (stats[,i]) <- "numeric"
  for (i in which (stats.mode == "boolean")) {mode (stats[,i]) <- "integer"; mode (stats[,i]) <- "logical"}
  ##
  stats <- as.list (stats)
  stats["filename"] <- file
  
  ##FEATURES
  features <- read.table (file = file, header = FALSE, sep = "\t", quote = "", colClasses = "character", nrows = 2, skip = 8)
  features <- features[,-1]
  ##
  features.mode  <- unlist (features[1,])
  features.names <- unlist (features[2,])
  ##
  fieldType <-  features.mode 
  names (fieldType) <- features.names
  ##
  fieldType <- gsub (pattern = "text",    replacement = "character", fieldType)
  fieldType <- gsub (pattern = "integer", replacement = "integer",   fieldType)
  fieldType <- gsub (pattern = "float",   replacement = "numeric",   fieldType)
  fieldType <- gsub (pattern = "boolean", replacement = "logical",   fieldType)
  ##
  fieldType <- as.list (fieldType)
  fieldType["filename"] <- file
  
  ###RETURN
  return (list (feparams = feparams, stats = stats, fieldType = fieldType))
}

## setwd ("/home/david/programas/mislibrerias/agilent/exdata/agilent")

## dir ()
## hs <- readAgilentHeaderS ()
## class (hs)
## length (hs)
## names (hs)
## hs[[1]]
## hs[1]
## class (hs[[1]])
## dim (hs[[1]])

## class (hs[[3]])

## unlist (hs[[3]][1,])

## as.vector (hs[[3]])
## as.character (hs[[3]])
## unlist (hs[[3]])


## rs <- readAgilentHeaderS ("agilent1.txt")
## class (rs)
## length (rs)
## names (rs)
## rs[[1]]
## class (rs[["feparams"]])
## class (rs[[3]])
## unlist (rs[[3]])


## dim (rs[["feparams"]])
## rownames (rs[["feparams"]])
## rownames (rs)
## rs[["feparams"]]
