##averageDuplicates.r
##2008-07-30 dmontaner@cipf.es
##2010-01-27 dmontaner@cipf.es
##2010-10-09 dmontaner@cipf.es
##2010-10-13 dmontaner@cipf.es
##2013-11-29 dmontaner@cipf.es


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
##' @param mat a numerical matrix.
##' @param ids character vector of length equal to the number of rows of \code{mat}.
##' @param original.rows logical indicating whether the mapping
##' between \code{ids} and original row position is to be kept.
##' @param summarystat the summary statistic: "mean" or "median"
##' @param na.rm should NA values be removed before computation?
##' @param verbose verbose
##' 
##' @return a matrix of averaged rows. \code{ids} are used as row names.
##' 
##' @examples
##' X <- data.frame(go=c("a","b","a","b"), exp1=c(3,4,5,7),
##'                 exp2=c(8,9,15,20), stringsAsFactors = FALSE)
##' X
##' M <- as.matrix (X[,c("exp1", "exp2")])
##' averageDuplicatedRows (mat = M, id = X[,"go"])
##' averageDuplicatedRows (mat = M, id = X[,"go"], original.rows = TRUE)
##' 
##' @export

averageDuplicatedRows <- function (mat, ids, original.rows = FALSE,
                                   summarystat = "mean", na.rm = TRUE, verbose = TRUE) {
 
  ##mat: matrix which rows are to be merged according to
  ##ids : vector indicating which rows are to be averaged
  ##summarystat: the summary statistic: "mean" or "median"
  
  if (nrow (mat) != length (ids)) stop ("nrow (mat) AND length (ids) DIFFER")
  if (!summarystat %in% c("mean", "median")) stop ("only \"mean\" or \"median\" are eligible summary statistics")
  
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
  
  ##dealing with unique ids rows (jut if there is more than one row)
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
    la1 <- lapply (1:length (tids), FUN = function (x) matrix (data = sp[[x]], nrow = tids[[x]]))
    if (summarystat == "mean") {
      la2 <- lapply (la1, FUN = colMeans, na.rm = na.rm)
    }
    if (summarystat == "median") {
      colMedians <- function (x, na.rm) {    ##NOT TOO EFFICIENT; it takes 10 times longer than colMeans
        apply (x, 2, median, na.rm = na.rm)
      }
      la2 <- lapply (la1, FUN = colMedians, na.rm = TRUE)
    }
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
    res <- list (means = res, originalrows = c(originalrows.u, originalrows.d))
  }
 
  ##output
  return (res)
}
