
###############################################################################  
##                                                                           ##  
## wNNSel - weighted nearest neighbor imputation using selected neighbors    ##  
##                                                                           ##  
## Author: Shahla Faisal shahla_ramzan@yahoo.com                             ##
##                                                                           ##   
###############################################################################  



##############
   artif.NAs <- function( x, testNA.prop, times.max  ) 
   {
  res <-  replicate(times.max, artifNA.cv(x, testNA.prop) , simplify=F)
  return(res)
            }

###################

cv.error.wNNSel = function(x, times.max, prelim.list, testNA.prop, k, lambda=0.1, method="2", c, m, kernel="gaussian")     # method=c("1","2")
 {
  if (!is.matrix(x)) stop("x should be a numeric data matrix")
  check(method, c, m)
  if(method=="1"){ tune=c }else{tune=m}

 MSIE = sapply(1:times.max, function(i, lambda, tune ) {
    x.train <-  prelim.list[[i]]$x.train
    remove.indices <- prelim.list[[i]]$remove.indices
    x <-  prelim.list[[i]]$x

    if(method=="2"){
    result <- wNNSel.impute(x.train, x.initial=x, lambda=lambda, method="2", m=tune, kernel=kernel, verbose=FALSE, verbose2=FALSE)  # m=m,
                   }

   if(method=="1"){
   result <- wNNSel.impute(x.train, x.initial=x, lambda=lambda, method="1", c=tune, kernel=kernel, verbose=FALSE, verbose2=FALSE)  # m=m,
                   }

    error = (result[remove.indices] - x[remove.indices])
    mean(error^2, na.rm=TRUE)
    }, lambda=lambda, tune=tune)

    MSE.final <-  mean(MSIE)
    if(method=="2"){
    result <- unlist(list( lambda=lambda,  m=tune, MSIE=MSE.final ) )  
    }

    if(method=="1") {
    result <- c( lambda=lambda,  c=tune, MSIE=MSE.final ) 
    }    
return(result)
}

###################
####    
impute.prelim = function(x, byrow = T, verbose=F) {
    missing.matrix = is.na(x)
    numMissing = sum(missing.matrix)
    if(verbose) {
      print(paste("imputing on", numMissing, "missing values with matrix size",
      nrow(x)*ncol(x), sep=" "))
    }
    if(numMissing == 0) {
      return ( list (missing.matrix = missing.matrix,
                numMissing = numMissing,
                missing.rows.indices = NULL,
                missing.cols.indices = NULL,
                x.missing = NULL) )
                           }
    missing.rows.indices = which(apply(missing.matrix, 1, function(i) {
          any(i)
        }))
    missing.cols.indices = which(apply(missing.matrix, 2, function(i) {
          any(i)
        }))
    if (byrow) x.missing = cbind(1:nrow(x),x)[missing.rows.indices,,drop=F]
    else x.missing = rbind(1:ncol(x),x)[,missing.cols.indices,drop=F]
    return ( list (missing.matrix = missing.matrix,
                   numMissing = numMissing,
                   missing.rows.indices = missing.rows.indices,
                   missing.cols.indices = missing.cols.indices,
                   x.missing = x.missing) )
 }

###################
#  1a function
cv.impute.prelim = function(x, test.fraction = 1/20) {
  n = nrow(x) * ncol(x)
  missing.matrix = is.na(x)
  valid.data = which(!missing.matrix)

  remove.indices = sample(valid.data, test.fraction*length(valid.data))
  x.train = x; x.train[remove.indices] = NA

  return (list(remove.indices = remove.indices,
               x.train = x.train))
}


###########

result.mat <- function(input.mat)
{
res.mat <- input.mat[,1:3]
for(  i in seq(1,ncol(input.mat),by=3)[-1]  )
{  res.mat <- rbind(res.mat, input.mat[,i:(i+2)])   }
return(res.mat)
}


###########

 impute.fn.wNNSel = function(values, distances, lambda, kernel, k, useAll=TRUE, verbosefn=FALSE)
        {
         if(useAll==TRUE)  {
                              k <- length(distances)
                              }else{
                                     if(!missing(k) & useAll==FALSE ) stopifnot( k <= length(distances) )
                                    }

          ranks = order(distances)
          smallest.distances = distances[ranks][1:k]
          smallest.distances = smallest.distances[!is.na(smallest.distances)]
          if ( is.numeric(smallest.distances) )   
          if (verbosefn) { print(" First k ordered smallest.distances are" );  print( smallest.distances) }
          knn.values = values[ranks][1:k]
          knn.values = knn.values[!is.na(knn.values)]
          if (verbosefn) { print("knn.values corresponding to smallest.distances are" ); print(knn.values) }

          if(is.null(kernel) )
          { knn.weights <- c( rep( 1/k, length(knn.values) ) )
             } else {
              knn.weights <- kernel.weight(smallest.distances, lambda, kernel )
                   }

          if (verbosefn) {  print("corresponding knn.weights are" ) ; print (knn.weights) }
          estimated.value <-  sum(knn.values * knn.weights)          
          if (verbosefn) { print("The imputed value is" ); print(estimated.value) }
          return(estimated.value)
             }                   ## end of impute.fn



###################

compute.d1c <- function(x, R) 
 {
  prelim = impute.prelim(x)
  if (prelim$numMissing == 0) { stop("x matrix has no missing values"); return (x)}
  dist.mat <- matrix(nrow=0,ncol=nrow(x)+2)
  colnames(dist.mat) <- c("row.index","col.index",rep("",nrow(x)) )

  for(i in 1:nrow(x))
                {
                if(  sum( is.na(x[i,]) ) > 0  )
                      {
                       miss.index <- which(is.na(x[i,]))

                       for(j in 1:length(miss.index))
                            {
                              temp<-matrix(nrow=nrow(x), data= rep( x[i,], nrow(x)),byrow=T )
                              dist.vec <- matrix(data= rowMeans(  abs(temp-x)* matrix( nrow=nrow(x),ncol=ncol(x),data=R[miss.index[j],],byrow=T) ,na.rm=T),byrow=T)
                              dist.mat <- rbind(dist.mat, c(i,miss.index[j],dist.vec) )
                            }       
                      }  
                }   
     return(dist.mat)
      }
###################
 compute.d2c  <-  function(x, R) {
    prelim = impute.prelim(x)
    if (prelim$numMissing == 0) { stop("x matrix has no missing values"); return (x)}
    dist.mat <- matrix(nrow=0,ncol=nrow(x)+2)
    colnames(dist.mat) <- c("row.index","col.index",rep("",nrow(x)) )
                for(i in 1:nrow(x))
                {
                if(  sum( is.na(x[i,]) ) > 0  )
                      {
                       miss.index <- which(is.na(x[i,]))

                       for(j in 1:length(miss.index))
                            {
                              temp<-matrix(nrow=nrow(x), data= rep( x[i,], nrow(x)),byrow=T )
                              dist.vec <- matrix(data= sqrt(rowMeans(  ( (temp-x)^2 )* matrix( nrow=nrow(x),ncol=ncol(x),data=R[miss.index[j],],byrow=T) ,na.rm=T) ),byrow=T)
                              dist.mat <- rbind(dist.mat, c(i,miss.index[j],dist.vec) )
                            } 
                      }   
                }  
     
          return(dist.mat)                            
         }
###################             
 compute.d1   <-  function(x) {
    prelim = impute.prelim(x)
    if (prelim$numMissing == 0) { stop("x matrix has no missing values"); return (x)}
    dist.mat <- matrix( nrow=nrow(x),ncol=nrow(x))
       for(i in 1:nrow(x) )
        {
         temp<-matrix(nrow=nrow(x), data= rep( x[i,], nrow(x)),byrow=T )
         dist.mat [i,]<- matrix(data= rowMeans(abs(temp-x),na.rm=T),byrow=T)
        }
     return(dist.mat)                            
    }
    
###################
 compute.d2   <-  function(x) {
    prelim = impute.prelim(x)
    if (prelim$numMissing == 0) { stop("x matrix has no missing values"); return (x)}            

     dist.mat <- matrix( nrow=nrow(x),ncol=nrow(x))
       for(i in 1:nrow(x) )
        {
          temp<-matrix(nrow=nrow(x), data= rep( x[i,], nrow(x)),byrow=T )
          dist.mat [i,]<- matrix( data=sqrt( rowMeans( (temp-x)^2, na.rm=T )  )   ,byrow=T )
        }
  return(dist.mat)                            
    }
   

###################

 check <- function( method, c, m )
        {
  if(method=="1" & missing(c) ) stop("c is required when method=1 ")
  if(method=="2" & missing(m) ) stop("m is required when method=2 ")
         }

###################






             