library(corrmi)
exprdata_set <- constructExprMatFromFile("/home/khaled/BioInformatics/GeneRegulatoryNetwork/MyTool/corrmi/data-raw/insilico_size10_1_multifactorial.tsv")
generank_mat <-
  matrix(-1,length(exprdata_set[1,]),length(exprdata_set[1,]))

l <- labels(exprdata_set)
rownames(generank_mat) <- l[[2]]
colnames(generank_mat) <- l[[2]]

genecorr_mat<-corrMat(exprdata_set)
#cmi_mat<-cmiMat(exprdata_set)
# for (i in 1:length(genecorr_mat[1,]))
# {
#   acol <- genecorr_mat[,i];
#   sortcol <- sort(acol,decreasing = TRUE)
#   generank_mat[,i] <- labels(sortcol)
# }
alpha<-0.1
lamda<-0.01

#for(i in 1:99){
#exprmat<-alpha * genecorr_mat + (1-alpha)*cmi_mat
#diag(exprmat)<- 0
#exprmat<-simp_Norm(exprmat)
#diag(exprmat)<- -1
GRN<-cmini(exprdata_set,lamda)
grn<-GRN$G
#   grn<-matrix(1,length(exprmat[1,]),length(exprmat[,1]))
#   for(j in 1:length(exprmat[1,])){
#     for(k in 1:length(exprmat[,1])){
#       if(exprmat[j,k]<lamda){
#         if(exprmat[j,k]>exprmat[k,j])
#         {
#           grn[j,k]<-0
#           exprmat[j,k]<-0
#         }else{
#           grn[k,j]<-0
#           exprmat[k,j]<-0
#         }
#       }
#     }
#   }
goldstd<-convertSortedRankTSVToAdjMatrix("/home/khaled/BioInformatics/GeneRegulatoryNetwork/MyTool/corrmi/data-raw/insilico_size10_1_goldstandard.tsv")
cnfsnMat<-evaluateGRN(grn,goldstd)
print(cnfsnMat)
cm<-matrix(c('TP','FP','FN','TN'),2,2)
print(cm)
tpr<-cnfsnMat[1,1]/(cnfsnMat[1,1]+cnfsnMat[1,2])
fpr<-cnfsnMat[2,1]/(cnfsnMat[2,1]+cnfsnMat[2,2])
print((cnfsnMat[1,1]+cnfsnMat[2,2])/sum(cnfsnMat))
alpha <- alpha + 1
