#' @title PackageLoaded
#' @description Determine if the package is loaded, if no package is loaded.
#' @param name  A character which is the name of package.
#' @return A boolean value.
#' @usage PackageLoaded(name)
#' @export


PackageLoaded<-function (name) {
  (paste("package:", name, sep = "") %in% search()) ||
    (name %in% loadedNamespaces())
}
