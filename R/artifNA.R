
###############################################################################
##                                                                           ##
## wNNSel - weighted nearest neighbor imputation using selected neighbors    ##
##                                                                           ##
## This R script contains the function to produce missing values in a given  ##
## data set completely at random.                                            ##
##                                                                           ##
## Author: Shahla Faisal shahla_ramzan@yahoo.com                             ##
##                                                                           ##
###############################################################################




#' Introduce MCAR Missing Values in a matrix 
#'
#' This function artificially introduces missing values in a data matrix under missing completely at random (MCAR) mechanism.
#'
#' @param x          a matrix, in which missing values are to be created.
#' @param miss.prop  proportion of missing values
#'
#' @return a matrix with missing values
#' @keywords NA
#' @export
#' @examples
#'  set.seed(3)
#'  x = matrix(rnorm(100),10,10)
#'  ## create 10% missing values in x
#'  artifNA(x, 0.10)


artifNA <- function(x, miss.prop=0.1 )
 {
    n <- nrow(x)
    p <- ncol(x)
    total <- n*p
    miss.ind <- sample(1:total, floor(total*miss.prop ) )
    x[miss.ind] <- NA
    return(x)
 }



