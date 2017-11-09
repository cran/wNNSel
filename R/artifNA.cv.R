

###############################################################################
##                                                                           ##
## wNNSel - weighted nearest neighbor imputation using selected neighbors    ##
##                                                                           ##
## This R script contains the function to produce missing values in a given  ##
## data set completely at random for cross validation.                       ##
##                                                                           ##
## Author: Shahla Faisal shahla_ramzan@yahoo.com                             ##
###############################################################################




#' Introduce MCAR Missing Values in a matrix for cross validation
#'
#' This function introduces additional missing values in a missing data matrix artificially.
#' The missing values are introduced under missing completely at random (MCAR) mechanism.
#' @param x          a matrix, in which missing values are to be created.
#' @param testNA.prop  proportion of missing values
#' @return  a list contatining a matrix with artifical missing values, removed indices and the provided x matrix
#' @seealso \code{\link{cv.wNNSel}}
#' @keywords NA cross-validation
#' @export
#' @examples
#'  set.seed(3)
#'  x = matrix(rnorm(100),10,10)
#'  ## create 10% missing values in x
#'  x.miss<- artifNA(x, 0.10)
#'  ## create another 10% missing values in x
#'  x.miss.cv<- artifNA.cv(x, 0.10)
#'  summary(x.miss)
#'  summary(x.miss.cv)


artifNA.cv <- function(x, testNA.prop=0.1 )
 {
    n <- nrow(x)
    p <- ncol(x)
    total <- n*p
  missing.matrix = is.na(x)
  valid.data = which(!missing.matrix)

  remove.indices = sample(valid.data, testNA.prop*length(valid.data))
  x.train = x
  x.train[remove.indices] = NA

  return (list(remove.indices = remove.indices, x.train = x.train, x=x))
}



