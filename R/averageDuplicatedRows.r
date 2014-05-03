##averageDuplicates.r
##2008-07-30 dmontaner@cipf.es
##2010-01-27 dmontaner@cipf.es
##2010-10-09 dmontaner@cipf.es
##2010-10-13 dmontaner@cipf.es
##2013-11-29 dmontaner@cipf.es
##2014-05-03 dmontaner@cipf.es

## min and max methods have not been tested much

##' @name averageDuplicatedRows
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @aliases insert
##' @keywords insert tiddler
##' 
##' @title Average duplicated rows
##'
##' @description Averages duplicated rows in a matrix.
##' 
##' @details The function averages rows in a matrix that are indexed by the same ids.
##' 
##' If original.rows is \code{TRUE} a length 2 list is returned.
##' The first element contains the row averaged values.
##' The second one keeps the link between ids and original row position.
##'
##' Some extra summary values are implemented apart from the "mean" value.
##' Available methods are "mean", "median", "min" and "max" for numerical matrices,
##' and "paste" and "pasteUnique" for character matrices.
##'
##' The "paste" and "pasteUnique" summaries use the function \code{paste}
##' for collapsing rows in character matrices.
##' "paste" collapses all values in the column (within each id) without handling duplicates.
##' "pasteUnique" removes duplicates and missing values before calling \code{paste}.
##'
##' @param mat a matrix. Generally a numerical one but see details.
##' @param ids character vector of length equal to the number of rows of \code{mat}.
##' @param original.rows logical indicating whether the mapping between \code{ids} and original row position is to be kept.
##' @param summarystat the summary statistic: "mean" or "median".
##' @param na.rm should NA values be removed before computation?
##' @param sep a character string to separate the terms when summarystat is "paste" or "pasteUnique".
##' @param verbose verbose
##' 
##' @return a matrix of averaged or collapsed rows. \code{ids} are used as row names.
##' 
##' @examples
##' X <- data.frame(go=c("a","b","a","b"), exp1=c(3,4,5,7),
##'                 exp2=c(8,9,15,20), stringsAsFactors = FALSE)
##' X
##' M <- as.matrix (X[,c("exp1", "exp2")])
##' by (M, X[,"go"], colMeans)
##' averageDuplicatedRows (mat = M, id = X[,"go"])
##' averageDuplicatedRows (mat = M, id = X[,"go"], original.rows = TRUE)
##'
##' T <- rbind (c("A", "aa", "aa", "aa"),
##'             c("B", "b1", "b2", "b3"),
##'             c("C", "cc", "cc", "cc"),
##'             c("B", "b1", "B2", "B3"),
##'             c("C", "cc", "cc", "cc"),
##'             c("C", "C1", NA, ""),
##'             c("D", "dd", "dd", "dd"))
##' T
##' averageDuplicatedRows (mat = T[,-1], id = T[,1], summarystat = "paste", sep = " * ")
##' averageDuplicatedRows (mat = T[,-1], id = T[,1], summarystat = "pasteUnique", sep = " * ")
##' 
##' @export

averageDuplicatedRows <- function (mat, ids, original.rows = FALSE,
                                   summarystat = "mean", na.rm = TRUE, sep = " /// ", verbose = TRUE) {
 
  ##mat: matrix which rows are to be merged according to
  ##ids : vector indicating which rows are to be averaged
  ##summarystat: the summary statistic: "mean", "median" ...

  mat <- as.matrix (mat)
    
  if (nrow (mat) != length (ids)) stop ("nrow (mat) AND length (ids) DIFFER")
  if (!summarystat %in% c("mean", "median", "min", "max", "paste", "pasteUnique")) stop ("only \"mean\", \"median\", \"min\", \"max\", \"paste\" or \"pasteUnique\" are eligible summary statistics")
  
  ##dealing with missing ids
  es.missing <- is.na (ids) | ids == ""
  conteo.missing <- sum (es.missing)
  if (conteo.missing > 0) {
    if (verbose) {
      warning (paste (conteo.missing, "missing ids where found in ids; they will all be tagged as \"missingID\""))
    }
    ids[es.missing] <- "missingID"
  }

  ##finding unique ids
  ids.count <- table (ids)[ids] #arranged as ids
  is.unique.ids <- ids.count == 1
  
  ##dealing with unique ids rows (just if there is more than one row)
  if (sum (is.unique.ids) > 1) {
    res.u <- mat[is.unique.ids,]
    rownames (res.u) <- ids[is.unique.ids]
  } else {
    res.u <- NULL
    is.unique.ids[is.unique.ids] <- FALSE
  }
  ##dealing with duplicated ids rows
  if (any (!is.unique.ids)) {
    sp  <- split (mat[!is.unique.ids,], ids[!is.unique.ids])
    tids <- as.list (table (ids[!is.unique.ids]))
    la1 <- lapply (1:length (tids), FUN = function (x) matrix (data = sp[[x]], nrow = tids[[x]]))  ## now we have a list of matrices spitted according to ids
    
    if (summarystat == "mean") {
      la2 <- lapply (la1, FUN = colMeans, na.rm = na.rm)
    }

    if (summarystat == "median") {
      colMedians <- function (x, na.rm) {    ##NOT TOO EFFICIENT; it takes 10 times longer than colMeans
        apply (x, 2, median, na.rm = na.rm)
      }
      la2 <- lapply (la1, FUN = colMedians, na.rm = na.rm)
    }

    if (summarystat == "min") {
      colMin <- function (x, na.rm) {
        apply (x, 2, min, na.rm = na.rm)
      }
      la2 <- lapply (la1, FUN = colMin, na.rm = na.rm)
    }
    
    if (summarystat == "max") {
      colMax <- function (x, na.rm) {
        apply (x, 2, max, na.rm = na.rm)
      }
      la2 <- lapply (la1, FUN = colMax, na.rm = na.rm)
    }
    
    if (summarystat == "paste") {
      colPaste <- function (x, sep) {
        apply (x, 2, paste, collapse = sep)
      }
      la2 <- lapply (la1, FUN = colPaste, sep = sep)
    }

    if (summarystat == "pasteUnique") {
      colPasteUnique <- function (x, sep) {
        pu <- function (x, sep = sep) {
          x <- x[!is.na (x)]
          x <- x[x != ""]
          x <- unique (x)
          x <- paste (x, collapse = sep)
          return (x)
        }
        apply (x, 2, pu, sep = sep)
      }
      la2 <- lapply (la1, FUN = colPasteUnique, sep = sep)
    }

    #names
    res.d <- matrix (unlist (la2), ncol = ncol (mat), byrow = TRUE)
    rownames (res.d) <- names (sp)
    colnames (res.d) <- colnames (mat)
  } else {
    res.d <- NULL
  }
  res <- rbind (res.u, res.d)
  
  ##recovering original positions
  if (original.rows) {
    if (sum (is.unique.ids) > 1) {
      originalrows.u <- which (is.unique.ids)
    } else {
      originalrows.u <- NULL      
    }
    if (any (!is.unique.ids)) {
      pos <- which (!is.unique.ids)
      pos.names <- names (pos)
      names (pos) <- NULL
      originalrows.d <- split (pos, pos.names)
    } else {
      originalrows.d <- NULL
    }
    res <- list (summarized = res, originalrows = c(originalrows.u, originalrows.d))
  }
 
  ##output
  return (res)
}
