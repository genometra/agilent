##readGPRHeaderS.r
##2009-08-23 dmontaner@cipf.es
##2010-10-24 dmontaner@cipf.es
##2013-11-29 dmontaner@cipf.es


##' @name readGPRHeaderS
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @keywords Read GPR header.
##' @seealso \code{\link{readAgilent}} \code{\link{readGPR}}
##' 
##' @references \url{http://www.moleculardevices.com/pages/software/gn_genepix_file_formats.html}
##' 
##' @title Read GPR headers.
##'
##' @description Reads header from GenePix GPR files.
##' 
##' @details Returns a list of two components 'params' and 'colnames'.
##' params is a dataframe having the header parameters; reads from files ares arranged in rows.
##' colnames is a character vector containing the names of the data columns.
##'
##' @param files names of the files which header is to be read.
##' @param samplenames names to be given to the arrays.
##' @param verbose verbose
##' @param includeFileNames should the filename be included?
##' 
##' @return A list of two components 'params' and 'colnames'.
##' 
##' @examples
##' setwd (file.path (system.file ("exampledata", package = "agilent")))
##' dir ()
##' 
##' hs <- readGPRHeaderS ()
##' class (hs)
##' length (hs)
##' hs$params[,1:5]
##' hs$params[,c("Wavelengths")]
##' hs$columns[1:5]
##'
##' hs <- readGPRHeaderS (files = c("GPR2.gpr", "GPR1.gpr"),
##'                       samplenames = c("second", "first"),
##'                       verbose = FALSE, includeFileNames = FALSE)
##' hs$param[,1:5]
##'
##' hs <- readGPRHeaderS ("GPR2.gpr")
##' class (hs)
##' hs$columns[1:5]
##'
## @import limma
##' @importFrom limma readGPRHeader
##' 
##' @export

readGPRHeaderS <- function (files = dir (pattern = ".gpr", ignore.case = TRUE),
                            samplenames = sub (".gpr", "", files, ignore.case = TRUE),
                            verbose = TRUE, includeFileNames = TRUE) {#2010-10-24 dmontaner@cipf.es

  ## INTERNAL FUNCTION #########################################################
  my.rbind <- function (r, s, stringsAsFactors = FALSE) {##2010-10-24 dmontaner@cipf.es
    ##appends DATAFRAMES with different column names
    ##not to be used with vectors or matrices
    if (is.null (r)) r <- as.data.frame (matrix (NA, nrow = 0, ncol = 0))
    if (is.null (s)) s <- as.data.frame (matrix (NA, nrow = 0, ncol = 0))
    ##
    r.names <- colnames (r)
    s.names <- colnames (s)
    i.names <- intersect (r.names, s.names)
    r.m.s <- setdiff (r.names, s.names)
    s.m.r <- setdiff (s.names, r.names)
    ##
    r <- cbind (r, matrix (NA, nrow = nrow (r), ncol = length (s.m.r), dimnames = list (rownames (r), s.m.r)), stringsAsFactors = stringsAsFactors)
    s <- cbind (s, matrix (NA, nrow = nrow (s), ncol = length (r.m.s), dimnames = list (rownames (s), r.m.s)), stringsAsFactors = stringsAsFactors)
    r <- rbind (r, s[,colnames (r)])
    ##
    return (r)
  }
  ##############################################################################
  
  ##require (limma)

  fileHeader <- readGPRHeader (files[1])
  nombres <- names (fileHeader)
  columnas <- read.table (file = files[1], header = FALSE, sep = "\t",
                          quote = "\"", colClasses = "character",
                          nrows = 1, skip = fileHeader$NHeaderRecords)

  headers <- NULL
  headers$params <- NULL
  ##headers$columns <- NULL
  for (fichero in files) {
    if (verbose) cat (paste ("reading header:", fichero), fill = TRUE)
    
    fileHeader <- readGPRHeader (fichero)
    fileColumns <- read.table (file = fichero, header = FALSE, sep = "\t",
                               quote = "\"", colClasses = "character",
                               nrows = 1, skip = fileHeader$NHeaderRecords)
    
    if (!identical (nombres, names (fileHeader))) warning (paste ("file", fichero, "seems to have a different gpr format. The header is not the same."))
    ##headers$params <- rbind (headers$params, as.data.frame (fileHeader, stringsAsFactors = FALSE)) ##falla aqui cuando los headers no son iguales
    headers$params <- my.rbind (headers$params, as.data.frame (fileHeader, stringsAsFactors = FALSE))

    if (!identical (columnas, fileColumns)) stop (paste ("file", fichero, "seems to have a different gpr format. Column names are not the same."))
  }
  
  headers$columns <- as.character (fileColumns)
  
  if (!is.null (samplenames)) {
    rownames (headers$params) <- samplenames
  }

  if (includeFileNames) {
    headers$params <- cbind (fileName = files, headers$params, stringsAsFactors = FALSE)
  }
  
  ##RETURN
  return (headers)
}
################################################################################

##readGPRHeader existe en limma
## readGPRHeader <- function (file) {#2009-07-20 dmontaner@cipf.es
  
##   ##reading the two first lines
##   dos <- readLines (con = file, n = 2)
##   dos <- gsub (" ", "", dos)
##   cuatro <- unlist (strsplit (dos, split = "\t"))
  
##   ##first element should be "ATF"
##   if (cuatro[1] != "ATF") stop (paste ("file", file, "does not start with ATF as GPR files do."))

##   ##creating header list
##   header <- list ()
##   ##header <- list (ATF = cuatro[2]) #version ##lo quitamos para que sea como limma

##   header[["NHeaderRecords"]] <- 2 + as.integer (cuatro[3]) #number of optional header records
##   header[["NDataFields"]] <- as.integer (cuatro[4]) #number of data columns
##   ##reading the remaining header lines
##   lineas <- readLines (con = file, n = header[["NHeaderRecords"]])
##   lineas <- lineas[-(1:2)]
##   lineas <- gsub ("\"", "", lineas)
##   ##
##   aux <- strsplit (lineas, split = "=")
##   info <- sapply (aux, function (x) x[2])  
##   names (info) <- sapply (aux, function (x) x[1])  
  
##   ##RETURN
##   header <- c(info, header)
##   return (header)
## }
