% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.wNNSel.R
\name{cv.wNNSel}
\alias{cv.wNNSel}
\title{Cross Validation for wNNSel Imputation}
\usage{
cv.wNNSel(x, kernel = "gaussian", x.dist = "euclidean", method = "2",
  m.values = seq(2, 8, by = 2), c.values = seq(0.1, 0.5, by = 0.1),
  lambda.values = seq(0, 0.6, by = 0.01)[-1], times.max = 5,
  testNA.prop = 0.05)
}
\arguments{
\item{x}{a \code{matrix} containing missing values}

\item{kernel}{kernel function to be used in nearest neighbors imputation. Default kernel function is "gaussian".}

\item{x.dist}{distance to compute, The default is \code{x.dist="euclidean"}
to compute Euclidean distance. Set \code{x.dist} to \code{NULL} to use Manhattan distance.}

\item{method}{convex function,  performs selection of variables. If \code{method="1"},
linear function is used and when if \code{method="c"}, power function is used.}

\item{m.values}{a \code{vector} of integer values, required when mehtod="2".}

\item{c.values}{a \code{vector} between 0 and less than 1. It is required when mehtod="1".}

\item{lambda.values}{a \code{vector}, for the tuning parameter \eqn{\lambda}}

\item{times.max}{maximum number of repititions for the cross validation procedure.}

\item{testNA.prop}{proportion of values to be deleted artificially for
cross validation in the missing matrix \code{x}. Default method uses 5 percent.}
}
\value{
a list containing 
   \item{lambda.opt}{optimal parameter selected by cross validation}
   \item{m.opt}{optimal parameter selected by cross validation}
   \item{MSIE.cv}{cross validation error}
}
\description{
This function aims to search for optimal values of the tuning parameters for the wNNSel imputation.
}
\details{
Some values are artificially deleted and wNNSel is run multiple times, varying  \eqn{\lambda} and  \eqn{m}.
For each pair of \eqn{\lambda} and \eqn{m}, compute MSIE on the subset of the data matrix x for which the
the values were deleted artificially. (See References for more detail).
}
\examples{
 set.seed(3)
 x.true = matrix(rnorm(100),10,10)
 ## create 10\% missing values in x
 x.miss = artifNA(x.true, 0.10)
 ## use cross validation to find optimal values
 result = cv.wNNSel(x.miss)
 ## optimal values are
 result$lambda.opt
 result$m.opt
 ## Now use these values to get final imputation
 x.impute = wNNSel.impute(x.miss, lambda=result$lambda.opt, m=result$m.opt)
 ## and final MSIE
 computeMSIE(x.miss, x.impute, x.true)
}
\author{
Shahla Faisal <shahla_ramzan@yahoo.com>
}
\references{
Tutz, G. and Ramzan,S. (2015). Improved methods for the imputation of missing data
by nearest neighbor methods.  \emph{Computational Statistics and Data Analysis}, Vol. 90, pp. 84-99.

Faisal, S. and Tutz, G. (2017). Missing value imputation for gene expression data by tailored nearest neighbors.
\emph{Statistical Application in Genetics and  Molecular Biology}.  Vol. 16(2), pp. 95-106.
}
\seealso{
\code{\link{artifNA.cv}}, \code{\link{wNNSel}}
}
\keyword{NA}
\keyword{cross-validation}
\keyword{wNNSel}
\keyword{weights}

