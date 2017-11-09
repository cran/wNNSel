
###############################################################################
##                                                                           ##
## wNNSel - weighted nearest neighbor imputation using selected neighbors    ##
##                                                                           ##
## This function converts the correlations to weights and thus performs      ##
##    selection of variables using a convex function.                        ##
##                                                                           ##
## Author: Shahla Faisal shahla_ramzan@yahoo.com                             ##
###############################################################################

convex.func <- function( r, method="2", m=2,  c=0.3)
  {
  
  if(method=="1" & missing(c) ) stop("c is required when method=1 ")
  if(method=="2" & missing(m) ) stop("m is required when method=2 ")
  
  if(method=="1"){
  convex <- ifelse(abs(r)<=c, 0, (abs(r)-c)/(1-c))
  }else{
        if(method=="2")    convex <- abs(r)^m
        }
  return(convex)
  }

