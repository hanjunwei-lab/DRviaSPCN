# DRviaSPCN

Drug Repurposing in Cancer via a Subpathway Crosstalk Network
> The `DRviaSPCN` is published in Bioinformatics. Please cite the following article when using `DRviaSPCN`: https://doi.org/10.1093/bioinformatics/btac611
        
        .

# Introduce

> We developed a novel software package (DRviaSPCN) that enables repurposing drugs via a subpathway (SP) crosstalk network. There are also several functions used to visualize the results such as visualization of the subpathway network structure of interest, chemical molecular formula of the drug or compound, and heatmap of the expression of subpathways in different sample types.

This network-based method consists of three major parts:

  - 1.Construction the SP network and calculation of the centrality scores of SPs to reflect the influence of SP crosstalk. 
  
  - 2.Calculating the enrichment scores of drug- and disease-induced dysfunctional SPs and weighted them by the centrality scores of SPs 
  
  - 3.Evaluating the drug-disease reverse association at the weighted SP level and identificating the cancer candidate drugs.


# How to install this package
> Installation method 1：
```
library(devtools); 
install_github("hanjunwei-lab/iATMEcell", build_vignettes = TRUE)
library(DRviaSPCN)
```
> Installation method 2：
```
install.packages("DRviaSPCN")
library(DRviaSPCN)
```


+  This package provides the `CalCentralityScore` function to calculate the centrality score of subpathways and the corresponding p-value.  

+  This package provides the `Optimaldrugs` function to calculate the DES of drugs and corresponding p-value.
  
+  This package provides the `plotSPW` function to plot SP network structure.

+  This package provides the `getMolecularFM` function to plot the chemical molecular formula of the drug or compound.

+  This package provides the `Disease2SPheatmap` function to plot heatmap of the activities of subpathways in different sample types that are regulated by disease.

+  This package provides the `Drug2SPheatmap` function to plot heatmap of the activities of subpathways in different sample types that are regulated by drugs.
  
+  This package provides the `GetExample` function to return example data and environment variables, such as gene expression profile, sample label and so on.</font>

> In addition, the essential data `DrugSPESCMatrix` and `DrugSPPvalueMatrix` which are subpathways weighted-ES induced by all drugs and statistic significance (p-value) of centrality score of subpathways regulated by all drugs were stored in our `DRviaSPCNData` package. Users could download and use this package by the following code:

```
library(devtools)
install_github("hanjunwei-lab/DRviaSPCNData",force = TRUE)
library(DRviaSPCNData)
### Get weighted-ES of subpathways
DrugSPESCMatrix<-GetData("DrugSPESCMatrix")
### Get p-value of subpathways centrality score
DrugSPPvalueMatrix<-GetData("DrugSPPvalueMatrix")
```

# Example 1 : Calculating the centrality scores of subpathways.
> The function `CalCentralityScore` was used to calculate the centrality scores of SPs to reflect the crosstalk influence, which were used as weights in the calculation of drug-disease reverse association score. 


```
###Load depend package
library(igraph)
###Obtain input data
GEP<-GetExample('GEP')# Get the gene expression profile
Slabel<-GetExample('Slabel')# Get the sample class label

###Run the function
CentralityScoreResult<-CalCentralityScore(ExpData=GEP,Label=Slabel,nperm = 1000)

```

# Example 2 : Calculating the drug-disease reverse association score and corresponding pvalue of drugs.
> The function `Optimaldrugs` is used to calculate the *DES* and statistic significance of drugs. The detailed algorithm can be seen in the introduction part. Users could screen out the optimal therapeutic drugs according to a specific threshold. Here we provide weighted and unweighted methods to calculate the score, which can be selected by parameters *weight = ''* . The screening criteria of the up- and down-regulated subpathways can be changed through the parameters *pcut = ''* and *topcut = ''*.

```
###Run the function
Opdrugresult<-Optimaldrugs(ExpData=GEP,Label=Slabel,DrugSPESC=DrugSPESCMatrix,
              CentralityScore=CentralityScoreResult,nperm=1000,topcut=10,
              pcut=0.01,weight=FALSE)
```


# Visualize 1: Plot a subpathway network structure graph.
> The function `plotSPW` used to plot a subpathway network structure graph. The user just needs to input an interest subpathway id such as "00020_4".

```
###load depend package
library(igraph)
###plot network graph of the subpathway "00020_4"
plotSPW("00020_4")
```

# Visualize 2: Plot a chemical molecular formula of the drug or compound .</font>

> The function `getMolecularFm` can obtain a chemical molecular formula of the drug or compound. Then users could visualize the molecular formula through function "plot".

```
###Load depend package
library(ChemmineR)
library(rvest)
###Obtain molecular formula and visualize it
Mole_formula<-getMolecularFm(drugname ="methotrexate")
plot(Mole_formula)
```


# Visualize 3: Plot a heatmap of the subpathways that are regulated by disease.</font>

> The function `Disease2SPheatmap` plots a heat map of the subpathways that are regulated by disease. The input is the result of function `CalCentralityScore`, disease gene expression profile and sample class in the expression profile. We map subpathways to the disease gene expression through ssgsea to get a subpathway abundance matrix. Then we visualize the matrix by heatmap. Users could change the threshold that is used to screen significant subpathways through the param *pcut*.

```
###Load depend package
library(GSVA)
library(pheatmap)

###Run the function
Disease2SPheatmap(CentralityScore=CentralityScoreResult,ExpData=GEP,Label=Slabel,pcut=0.05,
                   bk=c(-2,2),cluster.rows=FALSE,cluster.cols=FALSE,
                   show.rownames=TRUE,show.colnames=FALSE,
                   col=c("navy","firebrick3"),cell.width=NA,
                   cell.height=NA,scale="row",fontsize=7,
                   fontsize.row=9,fontsize.col=10)
```


# Visualize 4: Plot heatmaps of the subpathways that are regulated by drugs.</font>

> The function `Drug2SPheatmap` plots heatmaps of the subpathways that are regulated by drugs. The function input is a character which is drug name, disease gene expression profile and sample class in the expression profile. We map subpathways to the disease gene expression through ssgsea to get a subpathway abundance matrix. Then we visualize the matrix by heatmap. Users could change the threshold that is used to screen significant subpathways through the param *pcut*. The result of this function is a list including all heatmaps.

```
###Load depend package
library(GSVA)
library(pheatmap)
###Run the function
Drug2SPheatmap(drugname="methotrexate_HL60_6_8.8e-06",
              DrugSPPvalue=DrugSPPvalueMatrix,ExpData=GEP,
              Label=Slabel,pcut=0.05,bk=c(-2,2),cluster.rows=FALSE,
              cluster.cols=FALSE,show.rownames=TRUE,
              show.colnames=FALSE,col=c("navy","firebrick3"),
              cell.width=NA,cell.height=NA,scale="row",
              fontsize=6,fontsize.row=9,fontsize.col=10)


```

# Citation
These codes and data are intended for research use only.

``DRviaSPCN`` has already been published in *Bioinformatics*. If you use ``DRviaSPCN`` or these codes in your publication, please cite the paper:

Wu J, Li X, Wang Q, Han J. DRviaSPCN: a software package for drug repurposing in cancer via a subpathway crosstalk network. Bioinformatics. 2022;38(21):4975-4977. https://doi.10.1093/bioinformatics/btac611
        






