###############################################################################  
##                                                                           ##  
## wNNSel - weighted nearest neighbor imputation using selected neighbors    ##  
##                                                                           ##  
##    Compute weights for weighted nearest neighbor Imputation               ##  
##                                                                           ##
##                                                                           ##  
## Author: Shahla Faisal shahla_ramzan@yahoo.com                             ##
##                                                                           ##   
###############################################################################  

kernel.weight <- function(distances, lambda=0.3, kernel="gaussian")
{
 if(kernel=="triangular")
       {
       u <- distances/lambda                                   # ; print(u1)
   #    u <- u1/max(u1)                                     # ; print(u)
       Kernel <-  ifelse( abs(u)<=1, (1-abs(u)), 0 )       # ; print(Kernel)
 } else{
   if(kernel=="gaussian") {
         u <- distances/lambda                                   #  ; print(u)
         Kernel <-   (1/sqrt(2*pi))*exp((-1/2)*u^2)  }      #  ;  print(Kernel)
 }
    weight <- Kernel/ sum(Kernel)                            # ; print(weight)
    return(weight)
}
