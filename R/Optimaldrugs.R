getKS<-function(s,l) {
  if(class(s)=='character'&class(l)=='character'){
    V <- match(s,l)
  }
  if(class(s)=='numeric'&class(l)=='character'){
    V <- match(names(s),l)
  }
  if(class(s)=='character'&class(l)=='numeric'){
    V <- match(s,names(l))
  }
  if(class(s)=='numeric'&class(l)=='numeric'){
    V <- match(names(s),names(l))
  }
  V <- V[!is.na(V)]
  V <- sort(V)
  n <- length(l)
  t <- length(V)
  j <- 1:t
  a <- j/t - V/n
  a <- max(a)
  b <- V/n - (j - 1)/t
  b <- max(b)
  if (a > b) {
    ks = a
  }
  else {
    ks = -b
  }
  return(ks)
}

CalculateSES<-function (labels.list, correl.vector = NULL) {
  tag.indicator <- labels.list
  no.tag.indicator <- 1 - tag.indicator
  N <- length(labels.list)
  Nh <- length(labels.list[labels.list == 1])
  Nm <- N - Nh
  correl.vector <- abs(correl.vector)
  sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
  norm.tag <- 1/sum.correl.tag
  norm.no.tag <- 1/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag -
                  no.tag.indicator * norm.no.tag)
  max.ES <- max(RES)
  min.ES <- min(RES)
  ES <- signif(ifelse(max.ES > -min.ES, max.ES, min.ES), digits = 5)
  return(ES)
}


#'@title Identifying the optimal drugs
#'@description Function "Optimaldrugs" used to identify the optimal drugs for specific disease.
#'@param ExpData A gene expression profile of interest (rows are genes, columns are samples).
#'@param Label A character vector consist of "0" and "1" which represent sample class in gene expression profile. "0" means normal sample and "1" means disease sample.
#'@param DrugSPESC A matrix with n rows and m columns. n is the number of subpathways and m is the number of all drugs. The values in this matrix is weighted enrichmentscore of subpathways induced by each drug. The users could obtain this matrix from our example data.
#'@param CentralityScore The result of function "CalCentralityScore".
#'@param nperm Number of random permutations (default: 1000).
#'@param topcut The parameter "topcut" represents the number of selected SPs from the top or bottom of the ranked SP list. The topcut defaults to 10.
#'@param pcut The parameter "pcut" represents the threshold of statistical significance for screen SPs. The pcut defaults to 0.01.
#'@param weight A boolean value determines the method for calculating the drug-disease association score of the drug. "weight=FALSE"(default): Similar to "CMap" (Lamb et al., 2006), no weight is needed. "weight=TRUE": KS random walk statistic with individualized subpathway activity score as weight was used to calculate the drug-disease reverse association score.
#'@return A dataframe with four columns which are "Drug"(drug names),"DES"(drug enrichment score),"p-value"(statistical significance),"FDR"(adjusted statistical significance).
#'@importFrom stats p.adjust
#'@importFrom clusterProfiler GSEA
#'@usage Optimaldrugs(ExpData,Label,DrugSPESC,CentralityScore,nperm=1000,
#'                    topcut=10,pcut=0.01,weight=FALSE)
#'@export
#'@examples
#'##Obtain input data
#'#Weighted enrichmentscore of subpathways induced by each drug were stored
#'#in package "DRviaSPCNData". "DRviaSPCNData" has been uploaded to the
#'#github repository.Users can download and install through "install_github"
#'#function and set parameter url="hanjunwei-lab/DRviaSPCNData".
#'#After installing and loading package "DRviaSPCNData",
#'#users can use the following command to get the data.
#'#DrugSPESCMatrix<-GetData('DrugSPESCMatrix')
#'CentralityScoreResult<-GetExample("CentralityScoreResult")
#'GEP<-GetExample("GEP")
#'Slabel<-GetExample("Slabel")
#'\donttest{#Run the function
#'Opdrugresult<-Optimaldrugs(ExpData=GEP,Label=Slabel,DrugSPESC=DrugSPESCMatrix,
#' CentralityScore=CentralityScoreResult,nperm=1000,topcut=10,pcut=0.01,weight=FALSE)}


Optimaldrugs<-function(ExpData,Label,DrugSPESC,CentralityScore,nperm=1000,topcut=10,
                       pcut=0.01,weight=FALSE){
  havestats <- PackageLoaded("stats")
  if (havestats == FALSE) {
    stop("The 'stats' library, should be loaded first")
  }

  DE<-getDEscore(ExpData,Label)
  DE<-as.matrix(DE[which(abs(DE[,1])<100),])
  pathset<-data.frame(path='1',gene='1')
  for (i in 1:length(CentralityScore[,'SubPathID'])) {
    n<-CentralityScore[,'SubPathID'][i]
    g<-unlist(strsplit(CentralityScore[i,'Gene'],','))
    ps<-data.frame(path=rep(n,length(g)),gene=g)
    pathset<-rbind(pathset,ps)
  }
  pathset<-pathset[-1,]
  pathset[,1]<-as.character(pathset[,1])
  pathset[,2]<-as.character(pathset[,2])
  diseasegene<-DE[,1]
  names(diseasegene)<-rownames(DE)
  diseasegene<-sort(diseasegene,decreasing = T)


  diseaseES<-GSEA(diseasegene, TERM2GENE =pathset,exponent=1,minGSSize=0,
                  pAdjustMethod = "fdr", pvalueCutoff = 1)

  diseaseES<-diseaseES@result[,c(1,4,6)]

  diseasecentral<-as.numeric(CentralityScore[,'Centralscore'])
  names(diseasecentral)<-as.character(CentralityScore[,'SubPathID'])

  diseasecentral<-diseasecentral[match(rownames(diseaseES),names(diseasecentral))]
  a<-(sum(diseasecentral^2))^0.5
  b<-1+(diseasecentral/a)
  pathscore<-b*diseaseES[,2]
  result<-cbind(diseaseES[,1],pathscore)
  result1<-cbind(result,diseaseES[,3])
  result1<-as.data.frame(result1)
  result1[,1]<-as.character(result1[,1])
  result1[,2]<-as.numeric(as.character(result1[,2]))
  result1[,3]<-as.numeric(as.character(result1[,3]))
  colnames(result1)<-c('SubPathID','Weighted-ES','Pvalue')
  SubPathscore<-result1

  SubPathscore<-SubPathscore[order(SubPathscore[,2],decreasing = T),]

  uppath <- union(SubPathscore[SubPathscore[, 2] > 0 & SubPathscore[,3] < pcut, 1], SubPathscore[1:topcut,1])
  downpath <- union(SubPathscore[SubPathscore[, 2] < 0 & SubPathscore[,3] < pcut, 1],
                    SubPathscore[(length(SubPathscore[,1])-topcut+1):length(SubPathscore[,1]),1])
  if(weight==FALSE){
    ks_up <- apply(DrugSPESC, 2,function(y) {
      names(y) <- rownames(DrugSPESC)
      y <- sort(y,decreasing = TRUE)
      V <- match(uppath, names(y))
      V <- V[!is.na(V)]
      V <- sort(V)
      n <- length(y)
      t <- length(V)
      j <- 1:t
      a <- j/t - V/n
      a <- max(a)
      b <- V/n - (j - 1)/t
      b <- max(b)
      if (a > b) {
        ks_up = a
      }
      else {
        ks_up = -b
      }
      return(ks_up)
    })
    ks_down <- apply(DrugSPESC, 2,function(y) {
      names(y) <- rownames(DrugSPESC)
      y <- sort(y,decreasing = TRUE)
      V <- match(downpath, names(y))
      V <- V[!is.na(V)]
      V <- sort(V)
      n <- length(y)
      t <- length(V)
      j <- 1:t
      a <- j/t - V/n
      a <- max(a)
      b <- V/n - (j - 1)/t
      b <- max(b)
      if (a > b) {
        ks_down = a
      }
      else {
        ks_down = -b
      }
      return(ks_down)
    })
  }
  if(weight==TRUE){
    ks_up <- apply(DrugSPESC, 2,function(y) {
      names(y) <- rownames(DrugSPESC)
      y <- sort(y,decreasing = T)
      P.rank<-names(y)
      tag.indicator <- sign(match(P.rank,
                                  uppath, nomatch = 0))
      up.ES <- CalculateSES(tag.indicator, correl.vector = y)
      return(up.ES)
    })
    ks_down <- apply(DrugSPESC, 2,function(y) {
      names(y) <- rownames(DrugSPESC)
      y <- sort(y,decreasing = TRUE)
      P.rank<-names(y)
      tag.indicator <- sign(match(P.rank,
                                  downpath, nomatch = 0))
      down.ES <- CalculateSES(tag.indicator, correl.vector = y)
      return(down.ES)
    })
  }
  KS<-c()
  for (i in 1:length(ks_up)) {
    if(ks_up[i]*ks_down[i]>0){
      k<-0
    }else{
      k<-ks_up[i]-ks_down[i]
    }
    KS<-c(KS,k)
  }
  KS<-data.frame(KS,ks_up)

  #names(KS)<-colnames(DrugSPESC)
  ##标准化
  MA<-max(KS[,1])
  MI<-min(KS[,1])
  for(i in 1:length(KS[,1])){

    if(KS[,1][i]>0){
      KS[,1][i]<-KS[,1][i]/MA
    }
    if(KS[,1][i]<0){
      KS[,1][i]<--(KS[,1][i]/MI)
    }
  }
  KS<-KS[order(KS[,1],KS[,2],decreasing = TRUE),]
  KSorder<-KS[,1]
  names(KSorder)<-rownames(KS)
  ###按药名提取集合，富集到KS上，计算初始ks值
  drugname<-c()
  for (i in 1:length(KSorder)) {
    d<-unlist(strsplit(names(KSorder[i]),'_'))[1]
    drugname<-c(drugname,d)

  }
  uniquename<-unique(drugname)

  druglist<-list()
  for (i in 1:length(uniquename)) {
    druglist[[i]]<-KSorder[which(drugname==uniquename[i])]
    names(druglist)[i]<-uniquename[i]
  }

  ks_list <- lapply(druglist, getKS,KSorder)##计算初始ks值

  ##扰动
  permlist<-list()
  for (i in 1:length(druglist)) {
    permmatrix<-matrix(NA,length(druglist[[i]]),nperm)
    for (j in 1:nperm) {
      permmatrix[,j]<-sample(names(KSorder),length(druglist[[i]]),replace = FALSE)
    }
    permlist[[i]]<-permmatrix
  }

  pml<-lapply(permlist,function(x){
    apply(x, 2,getKS,KSorder)
  })##得到每个药扰动后的ks

  ##算p值
  pvalue<-c()
  for(i in 1:length(ks_list)){
    if(ks_list[[i]]>0){
      p<-length(which(pml[[i]]>ks_list[[i]]))/nperm
    }

    if(ks_list[[i]]<0){
      p<-length(which(pml[[i]]<ks_list[[i]]))/nperm
    }
    p<-length(which(pml[[i]]<ks_list[[i]]))/nperm
    pvalue<-c(pvalue,p)
  }
  fdr<-p.adjust(pvalue,'fdr',length(pvalue))
  result<-data.frame(matrix(unlist(ks_list),
                            nrow = length(ks_list),byrow = TRUE),
                     stringsAsFactors = FALSE)

  result<-cbind(names(ks_list),result)
  result<-cbind(result,pvalue)
  result<-cbind(result,fdr)
  result<-result[order(result[,3],decreasing = F),]
  colnames(result)<-c('Drug','DES','pvalue','FDR')
  result$DES<-round(result$DES,4)
  result$FDR<-round(result$FDR,3)
  return(result)
}
