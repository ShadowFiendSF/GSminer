% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/buildMapping.R
\docType{methods}
\name{CC2gene}
\alias{CC2gene}
\alias{CC2gene}
\alias{CC2gene,GeneGOMapping-method}
\title{Method CC2gene}
\usage{
\S4method{CC2gene}{GeneGOMapping}(object)
}
\arguments{
\item{object}{the GeneGOMapping object}
}
\value{
a list of GO to Gene mapping
}
\description{
Method CC2gene
}
\examples{

inputFile <- paste0(system.file(package = "GSminer"), "/extdata/TAIR.GO")
miner <- GSminer(inputFile = inputFile, sep = "\\t")
#for an example how to use:
#CC2gene(miner@ggmapping)

}
\author{
Li Zhaohong && Wu Zefeng
}
