
###############################################################################  
##                                                                           ##  
## wNNSel - weighted nearest neighbor imputation using selected neighbors    ##  
##                                                                           ##  
## Author: Shahla Faisal shahla_ramzan@yahoo.com                             ##
##                                                                           ##   
###############################################################################  


#' Weighted Nearest Neighbor Imputation of Missing Values using Selected Variables
#'
#' This package introduces new non-parametric tools for the imputation of missing values in high-dimensional data.   
#'   It includes weighted nearest neighbor
#'   imputation methods that use distances for selected covariates. The careful
#'   selection of distances that carry information about the missing values yields an imputation
#'   tool. It does not require pre-specified \eqn{k}, unlike other kNN methods. 
#'   It can be used to impute missing values in high-dimensional data when \eqn{n<p}. 
#' 
#' \tabular{ll}{ Package: \tab wNNSel\cr Version: \tab 0.1 \cr Date: \tab
#' 2017-11-08\cr Depends: \tab R (>= 2.10) \cr License: \tab GPL (>=
#' 2) }
#'
#' The main function of the package is \code{\link{wNNSel}} for implementing the nonparameteric procedure of nearest neighbors imputaiton.
#' See  \code{\link{wNNSel}} for more details.
#' 
#' 
#' @name wNNSel-package
#' @aliases wNNSel-package 
#' @docType package
#' @author Shahla Faisal  <shahla_ramzan@yahoo.com>   
#'
##### Maintainer: Shahla Faisal 

#' @references Tutz, G. and Ramzan,S*. (2015). Improved methods for the imputation of missing data
#' by nearest neighbor methods.  \emph{Computational Statistics and Data Analysis}, Vol. 90, pp. 84-99.
#'
#' Faisal, S.* and Tutz, G. (2017). Missing value imputation for gene expression data by tailored nearest neighbors.
#' \emph{Statistical Application in Genetics and  Molecular Biology}.  Vol. 16(2), pp. 95-106.
#'
#' @note *Author's Last name changed to \emph{Faisal} from \emph{Ramzan} in 2016.
#' @keywords package NA weights wNNSel
#' @import stats
NULL


