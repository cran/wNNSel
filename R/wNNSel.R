
###############################################################################  
##                                                                           ##  
## wNNSel - weighted nearest neighbor imputation using selected neighbors    ##  
##                                                                           ##  
## Author: Shahla Faisal shahla_ramzan@yahoo.com                             ##
##                                                                           ##   
###############################################################################  

# Weighted nearest neighbor imputation of missing values using selected variables
#' Imputatin using wNNSel method. 
#'
#' \code{'wNNSel'} is used to impute the missing values particularly in high dimensional data.
#' It uses a cross validation procedure for selecting the best values of the tuning parameters.
#' It also works when the samples are smaller than the covariates.  
#'
#'
#' For each sample, identify missinng features.  For each missing feature
#' find the nearest neighbors which have that feature.  Impute the missing
#' value using the imputation function on the \emph{selected} vector of values
#' found from the neighbors.
#' By default the \code{wNNSel} method automatically searches for optimal values for a given data matrix.
#'
#' The default method uses \code{x.dist="euclidean"} including selected covariates.
#' The specific distancs are computed using important covariates only.
#' If \code{mehtod="1"}, the linear function in absolute value of \eqn{r} is used, defined by
#' \deqn{\frac{|r|}{1-c} - \frac{c}{1-c},}
#' for \eqn{|r|>c}, and, 0 , otherwise. 
#' By default, the power function \eqn{|r|^m } is used when \code{mehtod="2"}. For more detailed discussion, see references.
#'
#'
#'
#'
#' @param x  a numeric data \code{matrix} containing missing values
#' @param k an optional, the number of nearest neighbors to use for imputation.
#' @param useAll \code{logical}. If \code{TRUE}, all \emph{available} neighbors are used for the imputation.
#' @param x.initial an optional. A complete data matrix e.g. using mean imputation of \code{x}. If provided, it will be used for the computation of correlations.
#' @param x.true  a matrix of true or complete data. If provided, \code{MSIE} will be returned in the results list.
#' @param x.dist distance to compute. The default is \code{x.dist="euclidean"}, that uses the Euclidean distance. Set \code{x.dist} to \code{NULL} for Manhattan distance.
#' @param kernel kernel function to be used in nearest neighbors imputation. Default kernel function is "gaussian".
#' @param impute.fn the imputation function to run on the length k vector of values for a missing feature. Defaults to a weighted mean of the neighboring values, weighted by the specified \code{kernel}. If not specified then wNN imputation will be used by default.
####### @param convex \code{logical}. whether selection of variables should be performed in computation of distance. The default is \code{TRUE}.
#' @param convex logical. If \code{TRUE}, selected variables are used for the computation of distance.  The default is \code{TRUE}.
#' @param c.values a \code{vector} between 0 and less than 1. It is required when mehtod="1".
#' @param m.values a \code{vector} of integer values, required when mehtod="2".
#' @param lambda.values  a \code{vector}, for the tuning parameter \eqn{\lambda}
#' @param times.max maximum number of repititions for the cross validation procedure.
#' @param testNA.prop  proportion of values to be deleted artificially for
#'        cross validation in the missing matrix \code{x}. Default method uses 5 percent.
#' @param method  convex function,  performs selection of variables. If \code{method="1"}, linear function is used and the power function is used when \code{method="2"}.
#' @param withinFolds \code{logical}. Use only if the neighbors/rows belong to particular folds/groups. Default is set to \code{FALSE}.
#' @param folds a \code{list} of vectors specifying folds/groups for neighbors. lenght of list is equal to the number of folds/groups. Each element/vector of the list indicates row indices belonging to that particular group/fold.
#' @param verbose  logical. If \code{TRUE}, prints status updates
#' @keywords wNNSel NA  cross-validation
#' @seealso \code{\link{cv.wNNSel}}, \code{\link{wNNSel.impute}}
#' @return a list containing imputed data matrix, and cross validation results
#'    \item{x.impute}{imputed data matrix}
#'    \item{MSIE}{True error. Note it is only available when x.true is provided.}
#'    \item{lambda.opt}{optimal parameter selected by cross validation}
#'    \item{m.opt}{optimal parameter selected by cross validation}
#'    \item{MSIE.cv}{cross validation error}
#'

#' @references Tutz, G. and Ramzan,S. (2015). Improved methods for the imputation of missing data
#' by nearest neighbor methods.  \emph{Computational Statistics and Data Analysis}, Vol. 90, pp. 84-99.
#'
#' Faisal, S. and Tutz, G. (2017). Missing value imputation for gene expression data by tailored nearest neighbors.
#' \emph{Statistical Application in Genetics and  Molecular Biology}.  Vol. 16(2), pp. 95-106.


#' @examples
#'  set.seed(3)
#'  x.true = matrix(rnorm(100),10,10)
#'  ## create 10% missing values in x
#'  x.miss = artifNA(x.true, 0.10)
#'  ## imputed matrix
#'  result <- wNNSel(x.miss)
#'  result$x.impute
#'  ## cross validation result can be accessed using
#'  result$cross.val
#  # lambda.opt   m.opt    MSIE.cv
#  # 0.0900000  8.0000000  0.5598143

#  ## another example when true x is known
#  result2 <- wNNSel(x.miss, x.true=x.true , method="1")
#  ## The true MSIE
#  result2$MSIE
###  [1] 1.034372

#' @export

#########    start of wNNSel function       ########


  wNNSel <-  function( x,  x.initial=NULL, x.true=NULL, k, useAll=TRUE, x.dist="euclidean", kernel="gaussian",  method="2" ,  impute.fn,
                       convex=TRUE, m.values = seq(2,8, by=2),   c.values = seq(.1,.5, by=0.1) , lambda.values = seq(0,.6,by=.01)[-1] ,
                       times.max = 5, testNA.prop = 0.05 , withinFolds=FALSE, folds, verbose=TRUE )
    {
    if (!is.matrix(x) || !is.numeric(x) ) stop("x should be a numeric data matrix")

    if( !missing(k) && k >= nrow(x))      stop("k must be less than the number of rows in x")       

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

     if (verbose) print("Cross validation in process...")
     result = cv.wNNSel(x, kernel, x.dist, method, m.values, c.values, lambda.values, times.max,  testNA.prop  )
     if (verbose) print("Cross validation complete")

if(method=="2"){
lambda.opt <- result$lambda.opt
     m.opt <- result$m.opt
  x.impute =  wNNSel.impute(x, k, useAll=TRUE, x.initial=NULL, x.dist="euclidean",  kernel="gaussian",  lambda=lambda.opt,  convex=TRUE, method="2", m=m.opt, verbose=FALSE, verbose2=FALSE)
   cv.res=list( lambda.opt=lambda.opt, m.opt=m.opt, MSIE.cv=result$MSIE.cv  )
             } else{

if(method=="1"){
lambda.opt = result$lambda.opt
     c.opt = result$c.opt
  x.impute =  wNNSel.impute(x, k, useAll=TRUE, x.initial=NULL, x.dist="euclidean",  kernel="gaussian",  lambda=lambda.opt,  convex=TRUE, method="1", c=c.opt, verbose=FALSE, verbose2=FALSE)
       cv.res=list(  MSIE.cv=result$MSIE.cv, lambda.opt, c.opt=c.opt  )
             }
           }
  if(!missing(x.true)) { if(!is.null(x.true))   MSIE <- computeMSIE(x, x.impute, x.true)} else{  MSIE=NULL }

  return( list( x.impute=x.impute, MSIE=MSIE, cross.val=unlist(cv.res)  ) )
    }



