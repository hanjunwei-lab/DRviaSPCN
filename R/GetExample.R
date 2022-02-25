#' @title Get example data
#' @description This function is used to achieve exxample data.
#' @param exampleData  A character, should be one of"GEP","Slabel",
#' "CentralityScoreResult","Opdrugresult" and "heatmap.list".
#' @return example data
#' @usage GetExample(exampleData)
#' @export


GetExample<-function(exampleData){
  if(!exists("envData")) {
    utils::data("envData",package="Disease2Drug")
  }

  if (exampleData=="GEP")
  {
    dataset<- get("GEP",envir=envData)
    return(dataset)
  }
  if (exampleData=="Slabel")
  {
    dataset<- get("Slabel",envir=envData)
    return(dataset)
  }
  if (exampleData=="SubPathwayInfo")
  {
    dataset<- get("SubPathwayInfo",envir=envData)
    return(dataset)
  }
  if (exampleData=="GoInfo")
  {
    dataset<- get("GoInfo",envir=envData)
    return(dataset)
  }
  if (exampleData=="Jaccardscore")
  {
    dataset<- get("Jaccardscore",envir=envData)
    return(dataset)
  }
  if (exampleData=="GoSubPconGene")
  {
    dataset<- get("GoSubPconGene",envir=envData)
    return(dataset)
  }
  if (exampleData=="SubPathwaymapdata")
  {
    dataset<- get("SubPathwaymapdata",envir=envData)
    return(dataset)
  }
  if (exampleData=="Drugs_CID")
  {
    dataset<- get("Drugs_CID",envir=envData)
    return(dataset)
  }
  if (exampleData=="CentralityScoreResult")
  {
    dataset<- get("CentralityScoreResult",envir=envData)
    return(dataset)
  }
  if (exampleData=="Opdrugresult")
  {
    dataset<- get("Opdrugresult",envir=envData)
    return(dataset)
  }
  if (exampleData=="heatmap.list")
  {
    dataset<- get("heatmap.list",envir=envData)
    return(dataset)
  }
}
