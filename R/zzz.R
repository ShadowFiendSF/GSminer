#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
    pkgVersion <- packageDescription(pkgname, fields="Version")
    msg <- paste0(pkgname, " v", pkgVersion, "  ", pkgname, "\n\n")

    citation <- paste0("If you use ", pkgname, " in published research, please cite our papers: \n\n", 
    	               "Zefeng Wu, Zhaohong Li and Ruolin Yang: GSminer: An R package for generating gold standard applied in gene functional network inference.2018 (in submission)\n\n")

    packageStartupMessage(paste0(msg, citation))
}

.onUnload <- function (libpath) {
  library.dynam.unload("GSminer", libpath)
}

