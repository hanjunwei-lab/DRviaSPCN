#'@title Calculating Log2 Fold Change of genes
#'@description Function "getDEscore" uses gene expression profile to calculate Log2 Fold Change of genes.
#'@param inexpData A gene expression profile of interest (rows are genes, columns are samples).The data in the expression profile is best not be log2 converted.
#'@param Label A character vector consist of "0" and "1" which represent sample class in gene expression profile. "0" means normal sample and "1" means disease sample.
#'@return A one-column matrix of Log2 Fold Change which rownames is gene.
#'@usage getDEscore(inexpData, Label)
#'@export

getDEscore <- function(inexpData, Label){
  expdata<-as.matrix(inexpData)
  test<-matrix(nrow=nrow(expdata),ncol=1)
  rownames(test)<-rownames(expdata)
  colnames(test)<-c("Log2FC")
  logfc.vector<-apply(expdata, 1, function(x){
    ind1 <- which(Label == 1)
    ind2 <- which(Label == 0)
    expdata1 <- x[ind1]
    expdata2 <- x[ind2]
    rmean1 <- mean(expdata1)
    rmean2 <- mean(expdata2)
    FC <- as.numeric(rmean1)/as.numeric(rmean2)
    logFC<-log2(FC)
    return(logFC)
  })
  test[,1]<-logfc.vector
  test <- as.matrix(test[!is.infinite(test[,1]),])
  test <- as.matrix(test[!is.na(test[,1]),])
  test <- as.matrix(test[!is.nan(test[,1]),])
  colnames(test)<-c("Log2FC")
  return(test)
}

