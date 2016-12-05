evaluateGRN<-function(expressionMat,goldMat){
  FN<-0
  FP<-0
  TP<-0
  TN<-0
  if((length(expressionMat[1,])==length(goldMat[1,])) && (length(expressionMat[,1])==length(goldMat[,1]))){
    print("Dimension OK!")
  }else{
    stop("Dimension Problem!")
  }

  for(i in 1:length(expressionMat[1,])){
    for(j in 1:length(expressionMat[,1]))
    {
      if(expressionMat[i,j]==1 && goldMat[i,j]==1)
      {
        TP<-TP + 1
      }else if(expressionMat[i,j]==1 && goldMat[i,j]==0){
        FP<-FP + 1
      }else if(expressionMat[i,j]==0 && goldMat[i,j]==1){
        FN<-FN + 1
      }else{
        TN<-TN + 1
      }
    }
  }
  #   p <- sum(goldMat)
  #   # The amount of negative links (these could be valid links and are not positive) [self regulating links are not allowed and are consired to be neither positive or negative]
  #   n <- (dim(goldMat)[1] * dim(goldMat)[2]) - p - sum(rownames(goldMat) %in% colnames(goldMat) )
  #   # The total amount of valid links (positive +negative)
  #   t <- p+n
  #   print(p)
  #   print(n)
  #   print(t)
  sprintf("TP=%d",TP)
  sprintf("FP=%d",FP)
  sprintf("FN=%d",FN)
  sprintf("TN=%d",TN)
  confusion<-c(TP,FP,FN,TN)
  confusionM<-matrix(confusion,2,2)
}

# calcAUROC <- function(fpr,tpr){
#
#   return(	trapz(fpr,tpr))
# }
#
# calcAUPR <- function(rec,prec,tpfn){
#
#   return(	trapz(rec,prec) / (1-1/tpfn))
# }
