
###############################################################################  
##                                                                           ##  
## wNNSel - weighted nearest neighbor imputation using selected neighbors    ##  
##                                                                           ##  
## Author: Shahla Faisal shahla_ramzan@yahoo.com                             ##
##                                                                           ##   
###############################################################################  

#'
#' Weighted Nearest Neighbor Imputation of Missing Values using Selected Variables
#'
##### Imputatin using wNNSel method. 
#'
#' This function imputes the missing values using user-spefied values of the tuning parameters. 
#' It also works when the samples are smaller than the covariates.  
#'
#' For each sample, identify missinng features.  For each missing feature
#' find the nearest neighbors which have that feature.  Impute the missing
#' value using the imputation function on the \emph{selected} vector of values
#' found from the neighbors.
#'
#'
#' @param x  a \code{matrix} containing missing values
#' @param k an optional, the number of nearest neighbors to use for imputation.
#' @param useAll \code{logical}. The default is \code{useALL=TRUE}, that is, all \emph{available} neighbors are used for the imputation.
#' @param x.initial an optional. A complete data matrix e.g. using mean imputation of \code{x}. If provided, it will be used for the computation of correlations.
#' @param x.dist distance to compute. The default is \code{x.dist="euclidean"}, that uses the Euclidean distance. Set \code{x.dist} to \code{NULL} for Manhattan distance.
#' @param kernel kernel function to be used in nearest neighbors imputation. Default kernel function is "gaussian".
#' @param lambda  \code{scaler}, a tuning parameter
#' @param impute.fn the imputation function to run on the length k vector of values for a missing feature. 
#' Defaults to a weighted mean of the neighboring values, weighted by the specified \code{kernel}. If not specified then wNN imputation will be used by default.
#' @param convex logical. If \code{TRUE}, selected variables are used for the computation of distance.  The default is \code{TRUE}.
#' @param m \code{scaler}, a tuning parameter required by the power function.
#' @param c \code{scaler}, a tuning parameter required by the linear function.  
#' @param method  convex function,  performs selection of variables. If \code{method="1"}, linear function is used and the power function is used when \code{method="2"}.
#' @param withinFolds \code{logical}. Use only if the neighbors/rows belong to particular folds/groups. Default is set to \code{FALSE}.
#' @param folds a \code{list} of vectors specifying folds/groups for neighbors. lenght of list is equal to the number of folds/groups. 
#'          Each element/vector of the list indicates row indices belonging to that particular group/fold.
#' @param verbose  logical. If \code{TRUE}, prints status updates
#' @param verbose2 logical. If \code{TRUE}, prints status updates with more detail
#' @keywords wNNSel NA weights
#' @seealso \code{\link{cv.wNNSel}},  \code{\link{wNNSel}}
#' @return imputed data matrix
#' @examples
#'   set.seed(3)
#'   x = matrix(rnorm(100),10,10)
#'   x.miss = x > 1
#'   x[x.miss] = NA
#'   wNNSel.impute(x)
#'   wNNSel.impute(x, lambda=0.5, m=2)
#' @export
#'
#########    start of wNNSel.impute function       ########

wNNSel.impute <- function(x, k, useAll=TRUE, x.initial=NULL, x.dist="euclidean",  kernel="gaussian",  lambda=0.3, impute.fn, convex=TRUE, method="2", m=2,c=0.3, 
                        withinFolds=FALSE, folds,  verbose=TRUE, verbose2=FALSE) {

    if (!is.matrix(x) || !is.numeric(x) ) stop("x should be a numeric data matrix")

    if( !missing(k) && k >= nrow(x))      stop("k must be less than the number of rows in x")        ## this is now unnecessary if count.NA is used

    if( withinFolds==TRUE && missing(folds) ) stop("The argument folds is missing")

    col.miss <- apply(x, 2, function(j) all(is.na(j)))
    row.miss <- apply(x, 1, function(i) all(is.na(i)))
  
    if (any(col.miss)) {
      cat("column(s)", which(col.miss), "are entirely missing.")
      stop("Please fix missing columns.")
    }
   
    if (any(row.miss)) {
    cat("row(s)", which(row.miss), "are entirely missing.")
    stop("Please fix missing rows.")   
     }

    check(method, c, m)
       
    prelim = impute.prelim(x)
    if (prelim$numMissing == 0) { print("x matrix has no missing values"); return (x)}
    missing.matrix = prelim$missing.matrix
    x.missing = prelim$x.missing
    missing.rows.indices = prelim$missing.rows.indices


 if (missing(impute.fn)) {
    impute.fn <- impute.fn.wNNSel        
  }

    if (verbose) print("Computing distance matrix...")
   
    dist.mat = dist.wNNSel( x=x, x.initial=x.initial, x.dist=x.dist, convex=convex, c=c, m=m, method=method )

    if (verbose) print("Distance matrix complete")

    x.missing.imputed = t(apply(x.missing, 1, function(i)
      {
        rowIndex = as.numeric(i[1])         
        i.original = unlist(i[-1])          
        if(verbose2) print(paste("Imputing row", rowIndex,sep=" "))
        missing.cols = which(missing.matrix[rowIndex,])  
        if(length(missing.cols) == ncol(x))
        warning( paste("Row",rowIndex,"is completely missing",sep=" ") )

              imputed.values = sapply(missing.cols, function(j)
                {
                    if(withinFolds)          
                      {
                      all.neighbor.indices = which(!missing.matrix[,j])
                      valid.ind = as.vector(na.omit (match(folds,  all.neighbor.indices )))
                      neighbor.indices <- all.neighbor.indices[-c(valid.ind)]
                      } else {
                      neighbor.indices = which(!missing.matrix[,j])
                            }
                  if (verbose2) { print(paste("Row",rowIndex,"and column", j, "is being imputed", sep=" ") )  }
                  if (verbose2) { print("All nearest neighbour indices having non-missing values in jth col are") ; print(neighbor.indices) }

                  if (!is.null(dist.mat))
                    {   
                        if(convex==TRUE)  {
                            knn.dist.all = dist.mat[dist.mat[,1] == rowIndex & dist.mat[,2] == j,  3:ncol(dist.mat)][neighbor.indices]
                            }   else {
                              knn.dist.all = dist.mat [ rowIndex , neighbor.indices ]            
                              }

                      knn.dist <-   knn.dist.all [ !is.nan(knn.dist.all)  ]               
                      if (verbose2) { print(" All 'valid' nearest neighbour distances are " ) ;  print(knn.dist) }
                     }
                   impute.fn(x[neighbor.indices,j],  knn.dist, lambda, kernel, k )
                 } )       

        i.original[missing.cols] = imputed.values      
        i.original
      } ) )     

    x[missing.rows.indices,] = x.missing.imputed
    x.missing.imputed <- c(x.missing.imputed)

       x.imputed = x
       return(x)
}                  ##end of wNNSel.impute function


