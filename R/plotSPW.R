
##' @title Plot subpathway network graph
##' @description The function plotSPW can visualize subpathway network graph.
##' @param subpathwayID Subpathway id .A character vector.
##' @param layout A matrix of x-y coordinates with two dims.
##' Determine the placement of the nodes for drawing a graph.
##' @param margin A numeric. The value is usually between -0.5 and 0.5, which is able to zoom in or out a subpathway graph. The default is 0.
##' @param vertex.label.cex A numeric vector of node label size.
##' @param vertex.label.font A numeric vector of label font.
##' @param vertex.label.dist A numeric vector of label dist.
##' @param vertex.size A numeric vector of Node size. See plot.igraph.
##' @param edge.arrow.width Edge arrow width. The default is 3.
##' @param edge.label.cex Edge label size.
##' @param vertex.label.color A vector of node label colors. The default is black.
##' @param vertex.color A vector of node colors. The default is the KEGG node color.
##' @param vertex.frame.color A vector of node frame color. The default is dimgray.
##' @param edge.color A vector of edge color. The default is dimgray
##' @param edge.label.color A vector of edge label color. The default is dimgray.
##' @param sub A character string of subtitle.
##' @param main A character string of main title.
##' @return a subpathway map
##' @export
##' @importFrom graphics plot
##' @importFrom igraph  plot.igraph
##' @examples
##' # load depend package
##' library(igraph)
##' # plot network graph of the subpathway "00020_4"
##' plotSPW("00020_4")




plotSPW<-function(subpathwayID,layout=NULL,margin = 0,
                         vertex.label.cex = 0.8, vertex.label.font = 0.8,vertex.label.dist=1,
                         vertex.size = 12,edge.arrow.width = 3, edge.label.cex =0.6,
                         vertex.label.color = "black",
                         vertex.color ="#F08080",vertex.frame.color = "dimgray",
                         edge.color = "grey70", edge.label.color = "dimgray",sub = NULL,
                         main = NULL){
  haveigraph <- PackageLoaded("igraph")
  if (haveigraph == FALSE) {
    stop("The 'igraph' library, should be loaded first")
  }
  subpathwaydata<-GetExample("SubPathwaymapdata")
  plot(subpathwaydata[[subpathwayID]],layout=layout,margin=margin,vertex.label.cex =vertex.label.cex ,vertex.label.font=vertex.label.font,
       vertex.size=vertex.size,vertex.label.dist=vertex.label.dist,edge.arrow.width=edge.arrow.width,
       edge.label.cex=edge.label.cex,vertex.label.color=vertex.label.color,vertex.color=vertex.color,
       vertex.frame.color=vertex.frame.color,edge.color=edge.color,edge.label.color=edge.label.color,sub=sub,main=main)
}
