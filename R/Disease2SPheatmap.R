##' @title Plot a heat map of the subpathways activity regulated by disease
##' @description The "Disease2SPheatmap" function plots a heat map of the subpathways that are regulated by disease.  We map subpathways to the disease gene expression through ssgsea to get a subpathway abundance matrix. Then we visualize the matrix by heatmap.
##' @param CentralityScore A dataframe with seven columns those are subpath ID, subpathway name, subpathway size, genes in subpathway, centralscore (eigenvector centrality), Pvalue and FDR.
##' @param ExpData A gene expression profile of interest.
##' @param Label A character vector consist of "0" and "1" which represent sample class in gene expression profile. "0" means normal sample and "1" means disease sample.
##' @param pcut A numeric value which represent threshold. Subpathways with p-value less than this threshold will be screened out and visualized.
##' @param bk A numeric vector that covers the range of values. Users could adjust color depth through this parameter.
##' @param cluster.rows Boolean values determining if rows should be clustered or hclust object.
##' @param cluster.cols Boolean values determining if columns should be clustered or hclust object.
##' @param show.rownames Boolean specifying if row names are be shown.
##' @param show.colnames Boolean specifying if column names are be shown.
##' @param col Vector of colors used in heatmap.
##' @param cell.width Individual cell width in points. If left as NA, then the values depend on the size of plotting window.
##' @param cell.height Individual cell height in points. If left as NA, then the values depend on the size of plotting window.
##' @param scale Character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are "row", "column" and "none".
##' @param fontsize Base fontsize for the plot (default: 10).
##' @param fontsize.row Fontsize for rownames (default: 10).
##' @param fontsize.col Fontsize for colnames (default: 10).
##' @return A heat map
##' @importFrom GSVA gsva
##' @importFrom pheatmap pheatmap
##' @importFrom grDevices colorRampPalette
##' @usage Disease2SPheatmap(CentralityScore,ExpData,Label,pcut=0.05,bk=c(-2,2),
##'                   cluster.rows=FALSE,cluster.cols=FALSE,show.rownames=TRUE,
##'                   show.colnames=FALSE,col=c("navy","firebrick3"),
##'                   cell.width=NA,cell.height=NA,scale="row",fontsize=7,
##'                   fontsize.row=9,fontsize.col=10)
##' @export
##' @examples
##' #Load depend package
##' library(GSVA)
##' library(pheatmap)
##' #Obtain input data (The "CentralityScoreResult" is the result of function "CalCentralityScore")
##' GEP<-GetExample('GEP')
##' Slabel<-GetExample('Slabel')
##' CentralityScoreResult<-GetExample('CentralityScoreResult')
##' #Run the function
##' Disease2SPheatmap(CentralityScore=CentralityScoreResult,ExpData=GEP,Label=Slabel,
##'                   pcut=0.05,bk=c(-2,2),cluster.rows=FALSE,cluster.cols=FALSE,
##'                   show.rownames=TRUE,show.colnames=FALSE,col=c("navy","firebrick3"),
##'                   cell.width=NA,cell.height=NA,scale="row",fontsize=7,
##'                   fontsize.row=9,fontsize.col=10)


Disease2SPheatmap<-function(CentralityScore,ExpData,Label,pcut=0.05,bk=c(-2,2),
                             cluster.rows=FALSE,cluster.cols=FALSE,show.rownames=TRUE,
                             show.colnames=FALSE,col=c("navy","firebrick3"),
                             cell.width=NA,cell.height=NA,scale="row",fontsize=7,
                             fontsize.row=9,fontsize.col=10){
  haveGSVA <- PackageLoaded("GSVA")
  havepheatmap <- PackageLoaded("pheatmap")
  if (haveGSVA == FALSE) {
    stop("The 'GSVA' library, should be loaded first")
  }
  if (havepheatmap == FALSE) {
    stop("The 'pheatmap' library, should be loaded first")
  }


  p<-CentralityScore[which(CentralityScore[,6]<pcut),]

  pl<-list()
  for (i in 1:length(p[,1])) {
    pl[[i]]<-strsplit(p[i,4],',')[[1]]
  }

  names(pl)<-p[,1]

  #library(GSVA)
  spw_matrix = gsva(as.matrix(ExpData), pl, method = "ssgsea",
                    kcdf = "Gaussian", abs.ranking = TRUE,min.sz=2)
  colnames(spw_matrix)[which(Label=='0')]='normal'
  colnames(spw_matrix)[which(Label=='1')]='disease'
  spw_matrix<-spw_matrix[,order(colnames(spw_matrix))]
  colann=data.frame(Sample=factor(rep(names(table(colnames(spw_matrix))),
                                      table(colnames(spw_matrix)))))
  samples<-paste(colnames(spw_matrix),1:length(spw_matrix[1,]))
  colnames(spw_matrix)<-samples
  rownames(colann)<-samples
  colnames(colann)<-'Group'

  bk1<-c(seq(bk[1],-0.1,by=0.1),seq(0,bk[2],by=0.1))
  #library(pheatmap)
  pheatmap(spw_matrix,
           scale = scale,
           breaks = bk1,
           color = colorRampPalette(c(col[1], "white", col[2]))(50),
           cluster_rows=cluster.rows,cluster_cols=cluster.cols,
           annotation_col =colann,
           show_rownames=show.rownames,show_colnames=show.colnames,
           cellwidth = cell.width,cellheight = cell.height,
           fontsize = fontsize,fontsize_row = fontsize.row,
           fontsize_col = fontsize.col,
           main = "Heatmap of the activities of significant subpathways regulated by the disease"

  )
}
