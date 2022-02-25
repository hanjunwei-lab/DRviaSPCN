random_network<-function(kegg_random){
  adj.final1<-as.matrix(kegg_random)
  graph1 = graph.adjacency(adj.final1,mode=c("undirected"), weighted=TRUE,add.rownames=T)
  temp1 = page.rank(graph1, vids=V(graph1), directed=FALSE, damping=0.90, weights=NULL)
  rank1 = temp1$vector
  rank2 = as.matrix(rank1)
  return(rank2)
}

go_p_score_row<-function(genecount,descore,path_size){
  score_row<-rep(0,path_size)
  for(j in 1:length(genecount)){
    if(genecount[j]!=""){
      gene<-unlist(strsplit(genecount[j], split = ","))
      location<-match(gene,descore[,1])
      de1<-descore[location,2]
      de_score1<-median(as.numeric(de1))
      if (!is.na(de_score1)) {
        score_row[j]<-de_score1
      }
    }
  }
  return(score_row)
}


##' @title Calculating eigenvector centrality of subpathways
##' @description The function "CalCentralityScore" is used to calculate the eigenvector centrality of subpathways.
##' @param ExpData A gene expression profile of interest (rows are genes, columns are samples).
##' @param Label A character vector consist of "0" and "1" which represent sample class in gene expression profile. "0" means normal sample and "1" means disease sample.
##' @param nperm Number of random permutations (default: 1000).
##' @return A dataframe with seven columns those are subpathway ID, subpathway name, subpathway size, genes in subpathway, centralscore (eigenvector centrality), Pvalue and FDR.
##' @importFrom igraph graph.adjacency
##' @importFrom igraph V
##' @importFrom igraph page.rank
##' @importFrom stats median
##' @usage CalCentralityScore(ExpData,Label,nperm=1000)
##' @export
##' @examples
##' library(igraph)
##' #Obtain input data
##' GEP<-GetExample('GEP')
##' Slabel<-GetExample('Slabel')
##' #Run the function
##' \donttest{CentralityScoreResult<-CalCentralityScore(ExpData=GEP,Label=Slabel,nperm=1000)}



CalCentralityScore <- function(ExpData,Label,nperm=1000){

  haveigraph <- PackageLoaded("igraph")
  havestats <- PackageLoaded("stats")
  if (haveigraph == FALSE) {
    stop("The 'igraph' library, should be loaded first")
  }
  if (havestats == FALSE) {
    stop("The 'stats' library, should be loaded first")
  }

  Subpathway<-GetExample('SubPathwayInfo')
  Go<-GetExample('GoInfo')
  Jaccard<-GetExample('Jaccardscore')
  Go_SubPath_gene<-GetExample('GoSubPconGene')

  DE<-getDEscore(ExpData,Label)
  DE<-as.matrix(DE[which(abs(DE[,1])<100),])
  kegg<-cbind(as.character(Subpathway[,"SubPathID"]),as.character(Subpathway[,"Gene"]))
  path_size<-length(kegg[,1])
  go<-as.matrix(Go)
  go_size<-length(go[,"Go_BP"])
  jaccard<-as.matrix(Jaccard)
  go_path_gene<-as.matrix(Go_SubPath_gene)
  DEscore<-cbind(names(DE[,1]),abs(DE[,1]))
  median_score<-matrix(0,nrow=go_size,ncol=path_size)
  for(k in 1:go_size){
    con_gene<-go_path_gene[k,]
    row<-go_p_score_row(con_gene,DEscore,path_size)
    median_score[k,]<-row
  }
  go_kegg<-median_score*jaccard
  colnames(go_kegg)<-kegg[,1]
  rownames(go_kegg)<-go[,1]
  #######构建subpath-subpath网
  edge<-as.matrix(go_kegg)
  edget<-t(edge)
  kegg_kegg<-edget%*%edge
  #######计算centra,
  adj.final<-as.matrix(kegg_kegg)
  graph = graph.adjacency(adj.final,mode=c("undirected"), weighted=TRUE,add.rownames=TRUE)
  temp = page.rank(graph, vids=V(graph), directed=FALSE, damping=0.90, weights=NULL)
  rank = temp$vector
  rank1 = as.matrix(rank)

  if (nperm == 0) {
    p <- rep(1, length(Subpathway[, 1]))
    fdr <- p
    allresult <- cbind(Subpathway, rank1)
    allresult <- cbind(allresult, p)
    allresult <- cbind(allresult, fdr)
    allresult[, 1] <- as.character(allresult[, 1])
    colnames(allresult) <- c("SubPathID", "SubPathway",
                             "Size", "Gene", "Centralscore",
                             "Pvalue", "FDR")
    result <- as.data.frame(allresult)
    rankresult <- result[order(result$SubPathID), ]
    rownames(rankresult) <- c(1:dim(kegg)[1])
  }else{
    iter<-nperm
    Centrality_Scores<-matrix(nrow=path_size,ncol=iter+1)
    real.centra<-rank1
    Centrality_Scores[,1]<-real.centra
    real.subname<-rownames(real.centra)
    subpathway<-as.character(colnames(kegg_kegg))
    for(i in 1:iter){
      per_subpathway<-sample(subpathway,replace = F)
      per_kegg_kegg<-kegg_kegg
      colnames(per_kegg_kegg)<-per_subpathway
      rownames(per_kegg_kegg)<-per_subpathway
      per_centra<-random_network(per_kegg_kegg)
      location<-match(real.subname,per_subpathway)
      per_centra1<-per_centra[location,]
      Centrality_Scores[,i+1]<-per_centra1
    }
    adj = as.matrix(Centrality_Scores)
    perm_rank = adj[,2:(iter+1)]
    perm_rank<-as.matrix(perm_rank)
    orig_rank = adj[,1]
    pval = matrix(data=NA, nrow=path_size, ncol=1)
    for ( j in 1:path_size ) {
      pval[j] = sum(perm_rank[j,] > orig_rank[j])/iter
    }
    p_padjust<-p.adjust(pval,method = "fdr")
    pa<-as.numeric(p_padjust)
    fdr<-round(pa,3)
    allresult<-cbind(Subpathway,rank1)
    allresult<-cbind(allresult,pval)
    allresult<-cbind(allresult,fdr)
    allresult[,1]<-as.character(allresult[,1])
    colnames(allresult)<-c("SubPathID","SubPathway","Size","Gene","Centralscore","Pvalue","FDR")
    result<-as.data.frame(allresult)
    rankresult<-result[order(result$`Pvalue`), ]
    rownames(rankresult)<-c(1:dim(kegg)[1])
  }
  rankresult$Centralscore<-round(rankresult$Centralscore,4)
  return(rankresult)
}
