
###############################################################################  
##                                                                           ##  
## wNNSel - weighted nearest neighbor imputation using selected neighbors    ##  
##                                                                           ##  
##  his function computes specific distances for the imputation of a         ##  
##  missing value using wNNSel method.                                       ##  
##                                                                           ##  
## Author: Shahla Faisal shahla_ramzan@yahoo.com                             ##
##                                                                           ##   
###############################################################################  


#### Distance function for wNNSel imputation
####
#### This function computes specific distances for the imputation of a missing value using \code{wNNSel} method.
####
#### The specific distancs are computed using important covariates only.
#### If \code{mehtod="1"}, the linear function in absolute value of \eqn{r} is used, defined by
#### \deqn{\frac{|r|}{1-c} - \frac{c}{1-c},}
#### for \eqn{|r|>c}, and, 0 , otherwise. 
#### By default, the power function \eqn{|r|^m } is used when \code{mehtod="2"}. For more detailed discussion, see references.
####
#### @param x a matrix containing missing values
#### @param x.initial an optional. A complete data matrix e.g. using mean imputation of \code{x}. If provided, it will be used for the computation of correlations.
#### @param x.dist distance to compute. The default is \code{x.dist="euclidean"}, that uses the Euclidean distance. Set \code{x.dist} to \code{NULL} for Manhattan distance.
#### @param convex logical. If \code{TRUE}, selected variables are used for the computation of distance.  The default is \code{TRUE}.
#### @param m \code{scaler}, a tuning parameter required by the power function.
#### @param c \code{scaler}, a tuning parameter required by the linear function.  
#### @param method  convex function,  performs selection of variables. If \code{method="1"}, linear function is used and the power function is used when \code{method="2"}.
#### @keywords distance weights
#### @references Tutz, G. and Ramzan,S. (2015). Improved methods for the imputation of missing data
#### by nearest neighbor methods.  \emph{Computational Statistics and Data Analysis}, Vol. 90, pp. 84-99.
####
#### Faisal, S. and Tutz, G. (2017). Missing value imputation for gene expression data by tailored nearest neighbors.
#### \emph{Statistical Application in Genetics and  Molecular Biology}.  Vol. 16(2), pp. 95-106.
####
#### @examples
####  set.seed(3)
####  x = matrix(rnorm(100),10,10)
####  x.miss = x > 1
####  x[x.miss] = NA
####  dist.wNNSel(x)
#### @export

dist.wNNSel <- function( x, x.initial=NULL, x.dist="euclidean", convex=TRUE, method="2",  m=2 , c=0.3 )
  {
   if(is.null(x.dist) )
     {
            if(convex==TRUE)           
              {
               check(method, c, m)
               if(is.null(x.initial))  {myx = x} else {myx=x.initial }
               R <- convex.func(r=cor(myx, use="pairwise.complete.obs",method="pearson"), method, m, c)
              dist.mat = compute.d1c(x, R)
              } else{       
               dist.mat = compute.d1(x)      
              }     
      } else { 

              if(x.dist=="euclidean")          
             {
             if (convex==TRUE)
             {
               check(method, c, m)
               if(is.null(x.initial))  {myx = x} else {myx=x.initial }
               R <- convex.func(r=cor(myx, use="pairwise.complete.obs",method="pearson"), method, m, c)
               dist.mat = compute.d2c(x, R)
             }

              else  {
                  dist.mat = compute.d2(x)                 
                  }
            }
    }          
  return(dist.mat)
 }

