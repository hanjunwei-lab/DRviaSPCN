% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Drug2SPheatmap.R
\name{Drug2SPheatmap}
\alias{Drug2SPheatmap}
\title{Plot a heat map of the subpathways activity regulated by drugs}
\usage{
Drug2SPheatmap(drugname="",DrugSPPvalue,ExpData,Label,pcut=0.05,
         bk=c(-2,2),cluster.rows=FALSE,cluster.cols=FALSE,show.rownames=TRUE,
         show.colnames=FALSE,col=c("navy","firebrick3"),
         cell.width=NA,cell.height=NA,scale="row",fontsize=7,
         fontsize.row=9,fontsize.col=10)
}
\arguments{
\item{drugname}{A character which represent interest drug name with specific concentration, cell line and duration.}

\item{DrugSPPvalue}{A matrix which colunms represent drugs and rows respresent subpathways. Values in this matrix is the pvalue of subpathways centrality score regulated by drugs.}

\item{ExpData}{A gene expression profile of interest.}

\item{Label}{A character vector consist of "0" and "1" which represent sample class in gene expression profile. "0" means normal sample and "1" means disease sample.}

\item{pcut}{A numeric value which represent threshold. Subpathways with p-value less than this threshold will be screened out and visualized.}

\item{bk}{A numeric vector that covers the range of values. Users could adjust color depth through this parameter.}

\item{cluster.rows}{Boolean values determining if rows should be clustered or hclust object.}

\item{cluster.cols}{Boolean values determining if columns should be clustered or hclust object.}

\item{show.rownames}{Boolean specifying if row names are be shown.}

\item{show.colnames}{Boolean specifying if column names are be shown.}

\item{col}{Vector of colors used in heatmap.}

\item{cell.width}{Individual cell width in points. If left as NA, then the values depend on the size of plotting window.}

\item{cell.height}{Individual cell height in points. If left as NA, then the values depend on the size of plotting window.}

\item{scale}{Character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are "row", "column" and "none".}

\item{fontsize}{Base fontsize for the plot (default: 10).}

\item{fontsize.row}{Fontsize for rownames (default: 10).}

\item{fontsize.col}{Fontsize for colnames (default: 10).}
}
\value{
A heat map
}
\description{
The function "Drug2SPheatmap" plots heatmaps of the subpathways that are regulated by drugs. We map subpathways to the disease gene expression through ssgsea to get a subpathway abundance matrix. Then we visualize the matrix by heatmap.
}
\examples{
#Load depend package
library(GSVA)
library(pheatmap)
##Obtain input data
#Statistic significance of centrality score of subpathways
#induced by each drug were stored in package "DRviaSPCNData".
#"DRviaSPCNData" has been uploaded to the github repository.
#Users can download and install through "install_github" function and
#set parameter url="hanjunwei-lab/DRviaSPCNData".
#After installing and loading package "DRviaSPCNData",
#users can use the following command to get the data:
#DrugSPPvalueMatrix<-GetData('DrugSPPvalueMatrix')
GEP<-GetExample('GEP')
Slabel<-GetExample('Slabel')

#Run the function
\donttest{Drug2SPheatmap(drugname = "methotrexate_HL60_6_8.8e-06",
                 DrugSPPvalue=DrugSPPvalueMatrix,
                 ExpData=GEP,Label=Slabel,pcut=0.05,bk=c(-2,2),
                 cluster.rows=FALSE,cluster.cols=FALSE,show.rownames=TRUE,
                 show.colnames=FALSE,col=c("navy","firebrick3"),
                 cell.width=NA,cell.height=NA,scale="row",fontsize=7,
                 fontsize.row=9,fontsize.col=10)}
}
