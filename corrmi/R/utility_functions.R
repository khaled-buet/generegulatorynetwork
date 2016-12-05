convertSortedRankTSVToAdjMatrix <- function(inputFilename=NULL,input=NULL,outputFilename=NULL){

  # print(getwd())
  # Read the table from file
  if(!is.null(inputFilename)){
    tbl <- read.table(inputFilename, header = TRUE)
  }
  else{
    tbl <- input
  }
  # Sort by second column
  tbl <- tbl[sort.list(tbl[,2]),]
  # Get the number of unique elements in the predictors
  pred <- unique(tbl[,1])
  # Get the number of unique elements in the targets
  targ <- unique(tbl[,2])
  # Pre allocate return matrix
  m <- matrix(0.0, length(pred),length(targ))
  # Set the col/row-names
  colnames(m) <- targ
  rownames(m) <- pred
  # Get the duplicates
  dups <- duplicated(tbl[,2])
  # Get the starIndices of another column
  startIndices <- which(FALSE == dups)
  for (i in 1:(length(startIndices)-1)){
    targetToAdd <- tbl[startIndices[i],2]
    predToAdd   <- tbl[startIndices[i]:(startIndices[i+1]-1),1]
    valuesToAdd <- tbl[startIndices[i]:(startIndices[i+1]-1),3]
    colIndex   <- which(colnames(m) %in% targetToAdd)

    tmp <- c()
    for (i in predToAdd){
      tmp <- c(tmp, which(rownames(m) == i))
    }
    rowIndexes <- tmp
    m[rowIndexes,colIndex] <- valuesToAdd


  }
  targetToAdd <- tbl[startIndices[length(startIndices)],2]
  predToAdd   <- tbl[startIndices[length(startIndices)]:length(tbl[,1]),1]
  valuesToAdd <- tbl[startIndices[length(startIndices)]:length(tbl[,1]),3]
  tmp <- c()
  for (i in predToAdd){
    tmp <- c(tmp, which(rownames(m) == i))
  }
  rowIndexes <- tmp
  colIndex   <- which(colnames(m) %in% targetToAdd)
  m[rowIndexes,colIndex] <- valuesToAdd

  if(!is.null(outputFilename)){
    write.table(m,outputFilename)
  }
  else{
    return(m)
  }
}
