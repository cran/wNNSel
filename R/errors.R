
###############################################################################  
##                                                                           ##  
## wNNSel - weighted nearest neighbor imputation using selected neighbors    ##  
##                                                                           ##  
##   Mean Squared Imputation Error                                           ##  
##    Mean Absolute Imputation Error                                         ## 
##    Normalized Root Mean Squared Imputatoin Error                          ##   
##                                                                           ##  
## Author: Shahla Faisal shahla_ramzan@yahoo.com                             ##
##                                                                           ##   
###############################################################################  


#' Mean Squared Imputation Error
#'
#' This function computes the mean squared imputation error for a given complete/true data matrix,
#' imputed data matrix and the data matrix with missing values.
#' @param x.miss     a \code{matrix},  having missing values
#' @param x.impute   an imputed data \code{matrix}. Note that it should not contain any missing values.
#' @param x.true     complete/true data \code{matrix}. Note that it should not contain any missing values.
#' @return  value of MSIE
#' @keywords error
#' @export
#' @examples
#'  set.seed(3)
#'  x.true = matrix(rnorm(100),10,10)
#'  ## create 10% missing values in x
#'  x.miss = artifNA(x.true, 0.10)
#'  ## impute using wNNSel method
#'  x.impute = wNNSel.impute(x.miss)
#'  computeMSIE(x.miss, x.impute, x.true)

computeMSIE <- function( x.miss, x.impute, x.true )
  {
    x.true=as.matrix(x.true)
    x.miss=as.matrix(x.miss)
    x.impute=as.matrix(x.impute)

    na.index = which(is.na(x.miss))
    mse    =  mean( ( x.true[na.index] - x.impute[na.index] )^2)
    return(mse)
 }




#' Mean Absolute Imputation Error
#'
#' This function computes the mean absolute imputation error for a given complete/true data matrix,
#' imputed data matrix and the data matrix with missing values.
#' @param x.miss     a \code{matrix},  having missing values
#' @param x.impute   an imputed data \code{matrix}. Note that it should not contain any missing values.
#' @param x.true     complete/true data \code{matrix}. Note that it should not contain any missing values.
#' @return  value of MSIE
#' @keywords error
#' @export
#' @examples
#'   set.seed(3)
#'   x.true = matrix(rnorm(100),10,10)
#'   ## create 10% missing values in x
#'   x.miss = artifNA(x.true, 0.10)
#'   ## impute using wNNSel method
#'   x.impute = wNNSel.impute(x.miss)
#'   computeMAIE(x.miss, x.impute, x.true)

#####   x.impute = wNNSel.impute(x.miss, lambda=0.5, m=2)

computeMAIE <- function( x.miss, x.impute, x.true)
  {
    x.true=as.matrix(x.true)
    x.miss=as.matrix(x.miss)
    x.impute=as.matrix(x.impute)

    na.index = which(is.na(x.miss))
    maie =  mean(abs( ( x.true[na.index] - x.impute[na.index] ) ))
    return(maie)
 }





#' Normalized Root Mean Squared Imputatoin Error
#'
#' This function computes the nrmalized root mean squared imputation error for a given complete/true data matrix,
#' imputed data matrix and the data matrix with missing values.
#' @param x.miss     a \code{matrix},  having missing values
#' @param x.impute   an imputed data \code{matrix}. Note that it should not contain any missing values.
#' @param x.true     complete/true data \code{matrix}. Note that it should not contain any missing values.
#' @return  value of MSIE
#' @keywords error
#' @export
#' @examples
#'   set.seed(3)
#'   x.true = matrix(rnorm(100),10,10)
#'   ## create 10% missing values in x
#'   x.miss = artifNA(x.true, 0.10)
#'   ## impute using wNNSel method
#'   x.impute = wNNSel.impute(x.miss)
#'   computeNRMSE(x.miss, x.impute, x.true)
#####   x.impute = wNNSel.impute(x.miss, lambda=0.5, m=2)



computeNRMSE <- function(  x.miss, x.impute, x.true )
  {
    x.true=as.matrix(x.true)
    x.miss=as.matrix(x.miss)
    x.impute=as.matrix(x.impute)
    na.index = which(is.na(x.miss))
    nrmse <- sqrt(mean((x.impute[na.index]-x.true[na.index])^{2})/var(x.true[na.index]))
    return(nrmse)
 }

