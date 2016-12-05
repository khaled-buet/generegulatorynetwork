corrMat <- function(exprdata_set) {

  genecorr_mat <-
    matrix(-1,length(exprdata_set[1,]),length(exprdata_set[1,]))

  l <- labels(exprdata_set)
  rownames(genecorr_mat) <- l[[2]]
  colnames(genecorr_mat) <- l[[2]]


  for (i in 1:length(exprdata_set[1,]))
  {
    x <- exprdata_set[,i]
    for (j in 1:length(exprdata_set[1,]))
    {
      if (i != j)
      {
        y <- exprdata_set[,j]
        cr <- cor(x,y,method = "pearson")
        genecorr_mat[i,j] <- cr
      }

    }
  }
  return(genecorr_mat)
}
