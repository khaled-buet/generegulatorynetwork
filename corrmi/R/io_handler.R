constructExprMatFromFile <-
  function(filename, rowFormat = FALSE, seperator = "\t", isLabeled = TRUE, labels =
             NULL,sampleNames = FALSE,normalize = TRUE) {
    DEFAULT.GENE.LABEL = "GENE_"
    m <-
      read.table(
        filename,sep = seperator,header = FALSE,as.is = TRUE,colClasses = "character"
      )

    if (is.null(labels)) {
      if (isLabeled) {
        if (rowFormat) {
          labels <- m[,1]
        }
        else		   {
          labels <- m[1,]
        }
      }
      else{
        labels <-
          paste(DEFAULT.GENE.LABEL,seq(from = 1,to = dim(m)[1]),sep = "")
      }
    }

    if (any(duplicated(colnames(m)))) {
      stop("Labels for genes are not unique")
    }

    if (rowFormat) {
      if (sampleNames) {
        m <- m[-1,]
      }
      if (isLabeled) 	{
        m <- m[,-1]
      }
    }
    else{
      if (sampleNames) {
        m <- m[,-1]
      }
      if (isLabeled) 	{
        m <- m[-1,]
      }
    }

    m <- as.matrix(m)
    m <- apply(m,2, as.numeric)

    if (rowFormat) {
      m <- t(m)
    }

    colnames(m) <- labels
    if (normalize) {
      m <- scale(m)
    }

    return(m)
  }

simp_Norm<-function (expDat){
  expDat<-expDat - min(expDat)
  ans<-matrix(0, nrow=nrow(expDat), ncol=ncol(expDat))
  for(i in seq(ncol(expDat))){
    ans[,i]<-expDat[,i]/sum(expDat[,i])
  }
  colnames(ans)<-colnames(expDat)
  rownames(ans)<-rownames(expDat)
  return(ans)
}
