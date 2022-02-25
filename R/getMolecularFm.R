#' @title Plot chemical molecular formula  of drugs
#' @description The function "getMolecularFm" outputs the chemical molecular formula  of a drug or compound . The results can be visualized by the "plot" function.
#' @param drugname A character string of drug name.
#' @param main An overall title for the chemical structure graph.
#' @param sub A sub title for the chemical structure graph.
#' @return Chemical molecular formula  of the drug or compound.
#' @importFrom rvest html_text
#' @importFrom ChemmineR read.SDFset
#' @importFrom xml2 read_html
#' @usage getMolecularFm(drugname = "", main = "", sub = "")
#' @export
#' @examples
#' ##Load depend package
#' library(ChemmineR)
#' library(rvest)
#' # Obtain molecular formula and visualize it.
#' Mole_formula<-getMolecularFm(drugname ="methotrexate")
#' plot(Mole_formula)



getMolecularFm<-function(drugname = "", main = "", sub = "") {
  haveChemmineR <- PackageLoaded("ChemmineR")
  havervest <- PackageLoaded("rvest")
  if (haveChemmineR == FALSE) {
    stop("The 'ChemmineR' library, should be loaded first")
  }
  if (havervest == FALSE) {
    stop("The 'rvest' library, should be loaded first")
  }
  drugname <- unlist(strsplit(drugname, "\\("))[1]
  drugname1 <- tolower(drugname)
  Drugs_CID <- GetExample("Drugs_CID")
  drugCid <- Drugs_CID[which(Drugs_CID[, 1] == drugname1),
                       2]
  drug_url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/",
                    drugCid, "/record/SDF/?record_type=2d&response_type=display",
                    sep = "")
  cw <- try(read_html(drug_url))
  if ("try-error" %in% class(cw)) {
    stop("Please ensure smooth network connection")
  }
  drugnr <- html_text(cw)
  drugnr <- strsplit(drugnr, "\n")
  drugnr <- unlist(drugnr)
  sdfset <- read.SDFset(drugnr)
  if (main == "") {
    sdfset@ID <- drugname
  }
  else {
    sdfset@ID <- main
  }
  return(sdfset)

}





