
###############################################################################  
##                                                                           ##  
## wNNSel - weighted nearest neighbor imputation using selected neighbors    ##  
##                                                                           ##  
##  This function aims to autimatically search for optimal values of the     ##  
##  tuning parameters that yield smallest MSIE for a missing data matrix.    ##  
##                                                                           ##  
## Author: Shahla Faisal shahla_ramzan@yahoo.com                             ##
##                                                                           ##   
###############################################################################  



#' Cross Validation for wNNSel Imputation
#'
#' This function aims to search for optimal values of the tuning parameters for the wNNSel imputation.
##### that yield smallest MSIE for a missing data matrix.
#'
#'
#' Some values are artificially deleted and wNNSel is run multiple times, varying  \eqn{\lambda} and  \eqn{m}.
#' For each pair of \eqn{\lambda} and \eqn{m}, compute MSIE on the subset of the data matrix x for which the
#' the values were deleted artificially. (See References for more detail).
##### (For details, Tutz and Ramzan(2015), Faisal and Tutz (2017) ).
#'
#' @param x  a \code{matrix} containing missing values
#' @param kernel kernel function to be used in nearest neighbors imputation. Default kernel function is "gaussian".
#' @param x.dist distance to compute, The default is \code{x.dist="euclidean"}
#'               to compute Euclidean distance. Set \code{x.dist} to \code{NULL} to use Manhattan distance.
#' @param method  convex function,  performs selection of variables. If \code{method="1"},
#'                linear function is used and when if \code{method="c"}, power function is used.
#' @param c.values a \code{vector} between 0 and less than 1. It is required when mehtod="1".
#' @param m.values a \code{vector} of integer values, required when mehtod="2".
#' @param lambda.values  a \code{vector}, for the tuning parameter \eqn{\lambda}
#' @param times.max maximum number of repititions for the cross validation procedure.
#' @param testNA.prop  proportion of values to be deleted artificially for
#'        cross validation in the missing matrix \code{x}. Default method uses 5 percent.
#' @return a list containing 
#'    \item{lambda.opt}{optimal parameter selected by cross validation}
#'    \item{m.opt}{optimal parameter selected by cross validation}
#'    \item{MSIE.cv}{cross validation error}
#' 
#' 
#' @author Shahla Faisal <shahla_ramzan@yahoo.com>
#' @references Tutz, G. and Ramzan,S. (2015). Improved methods for the imputation of missing data
#' by nearest neighbor methods.  \emph{Computational Statistics and Data Analysis}, Vol. 90, pp. 84-99.
#'
#' Faisal, S. and Tutz, G. (2017). Missing value imputation for gene expression data by tailored nearest neighbors.
#' \emph{Statistical Application in Genetics and  Molecular Biology}.  Vol. 16(2), pp. 95-106.
#' @seealso \code{\link{artifNA.cv}}, \code{\link{wNNSel}}
#' @keywords wNNSel NA weights cross-validation 
#' @export
#' @examples
#'  set.seed(3)
#'  x.true = matrix(rnorm(100),10,10)
#'  ## create 10% missing values in x
#'  x.miss = artifNA(x.true, 0.10)
#'  ## use cross validation to find optimal values
#'  result = cv.wNNSel(x.miss)
#'  ## optimal values are
#'  result$lambda.opt
#'  result$m.opt
#'  ## Now use these values to get final imputation
#'  x.impute = wNNSel.impute(x.miss, lambda=result$lambda.opt, m=result$m.opt)
#'  ## and final MSIE
#'  computeMSIE(x.miss, x.impute, x.true)



  
cv.wNNSel <-  function( x,   kernel="gaussian", x.dist="euclidean", method="2" ,
          m.values = seq(2,8, by=2),   c.values = seq(.1,.5, by=0.1) , lambda.values = seq(0,.6,by=.01)[-1] ,
         times.max = 5, testNA.prop = 0.05 )
 {

if(method=="1")   tune.vec = c.values
if(method=="2")   tune.vec = m.values

prelim.list = artif.NAs( x, testNA.prop , times.max )
               
mse.mat <-t( sapply(lambda.values, function(ii)
 {
 sapply(tune.vec, function (jj) {
 #res <- cv.error.wNNSel( testNA.prop=testNA.prop, k=k, lambda=ii, method="2", m=jj)
 if(method=="2") {  res <- cv.error.wNNSel(x, times.max=times.max, lambda=ii, method="2", m=jj, kernel=kernel, prelim.list=prelim.list)
     } else{  if(method=="1")   res <- cv.error.wNNSel(x, times.max=times.max, lambda=ii, method="1", c=jj, kernel=kernel , prelim.list=prelim.list)  }

 })
})    )

   lambda.tune.mse.mat <- result.mat(mse.mat)


lambda.opt <- lambda.tune.mse.mat[  which.min(lambda.tune.mse.mat[,3])  , 1 ]
tune.opt   <- lambda.tune.mse.mat[  which.min(lambda.tune.mse.mat[,3])  , 2 ]
mse.opt    <- lambda.tune.mse.mat[  which.min(lambda.tune.mse.mat[,3])  , 3 ]

 if(method=="1")  result <- list(lambda.opt=lambda.opt, c.opt=tune.opt, MSIE.cv=mse.opt)
 if(method=="2")  result <- list(lambda.opt=lambda.opt, m.opt=tune.opt, MSIE.cv=mse.opt)

return( result  )
 }

