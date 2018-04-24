#' Class "GeneGOMapping"
#' A S4 class for mapping between gene and GO terms
#' This class is the solt of the GSminer which used to select genes.
#' 
#' @name GeneGOMapping-class
#' @aliases GeneGOMapping-class 
#' @docType class
#' @slot inputFile the name of input file which contain the mapping 
#'       relationship between the gene and GO terms
#' @slot sep used in the input file by default is tab
#' @slot mapping is a private env
#' @export GeneGOMapping
#' @exportClass GeneGOMapping
#' @import methods
#' @examples
#'
#' inputFile <- paste0(system.file(package = "GSminer"), "/extdata/TAIR.GO")
#' ggmap <- GeneGOMapping(inputFile = inputFile, sep = "\t")
#' ###	
#' ggmap <- methods::new("GeneGOMapping", inputFile = inputFile, 
#'                         sep = "\t")	
#' 
#' @author Li Zhaohong && Wu Zefeng
#' @seealso \code{\linkS4class{GSminer}}
#' @keywords classes


GeneGOMapping <- setClass(Class = "GeneGOMapping",
						  slots = c(inputFile = "character", sep = "character", mapping = "environment"),
						  prototype = prototype(inputFile = inputFile <- paste0(system.file(package = "GSminer"), "/extdata/TAIR.GO"),
						  sep = "\t", mapping = new.env(parent=emptyenv()))
						  )

#
# initialize method for the GeneGOMapping class
# mimic the RAII conception in C++
#
setMethod("initialize", "GeneGOMapping", function(.Object, inputFile, sep, ...)
		  {
		        
			.Object@inputFile <- ifelse(missing(inputFile), inputFile <- paste0(system.file(package = "GSminer"), "/extdata/TAIR.GO"), inputFile)
			.Object@sep <- ifelse(missing(sep), "\t", sep)
			mapping <- new.env(parent=emptyenv())
			mapping$mapping <- list()
			.Object <- methods::callNextMethod(.Object, ...)
			return(.Object)
		  }
)

#'
#' Validation of GeneGOMapping object
#' @param object instance of GeneGOMapping for gene to GO terms mapping
#' @return TRUE for validation
#'

MappingFileValid <- function(object)
{
	if(length(object@inputFile)==1 && length(object@sep)==1)
	{
		if(!file.exists(object@inputFile))
			return(paste("Gene to GO mapping file", object@inputFile, "must exist!", sep=" "))
		if(identical(file.size(object@inputFile), 0))
			return(paste("Gene to GO mapping file", object@inputFile, "size must not be empty!", sep=" "))
		if(!TestSep(object@inputFile, object@sep))
			return(paste("Gene to GO mapping file", object@inputFile, "must hava correct separator!", sep=" "))
		TRUE
	}else{
		return(paste("insufficient parameters for the constructor!"))
	}
}

TestSep <- function(filename, sep)
{
	if(file.exists(filename))
	{
		con<-file(filename, open="r")
		while(length(oneLine<-readLines(con, n=1))>0)
		{
			if(!grepl("^#", oneLine, perl=TRUE))
			{
				if(length(unlist(strsplit(oneLine, split=sep, perl=TRUE)))>=2)
				{
					close(con)
					return(TRUE) 
				}else{
					close(con)
					return(FALSE)
				}
			}
		}
		close(con)
	}
	FALSE
}

#
# Validation of GeneGOMapping object
#
setValidity("GeneGOMapping", MappingFileValid)

#' Class "GSminer"
#' 
#' GSminer has three slots:
#' 1. ggmapping object for gene to GO terms mapping
#' 2. BPoffsprings a list for BP GO terms to gene mapping
#' 3. CCoffsprings a list for CC GO terms to gene mapping
#' Note: the BPoffspring and CCoffspring mimic the behavior of pass by const reference in C++
#' @name GSminer-class
#' @aliases GSminerClass
#' @docType class
#' @slot ggmapping which is the instance of GeneGOMapping for gene to GO terms mapping
#' @slot mapping is a private env
#' @param inputFile input file name(with absolute path)
#' @param sep input file separator by default is tab
#' @export GSminer
#' @exportClass GSminer
#' @import methods
#' @examples
#'
#' inputFile <- paste0(system.file(package = "GSminer"), "/extdata/TAIR.GO")
#' miner <- GSminer(inputFile = inputFile, sep = "\t")
#' ggmap <- methods::new("GeneGOMapping", inputFile = inputFile, sep = "\t")
#' miner@ggmapping <- ggmap
#' loadBPoffsprings(miner)
#'
#' @seealso \code{\linkS4class{GeneGOMapping}}
#' @author Li Zhaohong && Wu Zefeng
#' @keywords classes


GSminer <- setClass(Class = "GSminer",
		              slots = c(ggmapping = "GeneGOMapping", mapping = "environment"),
					  prototype = prototype(ggmapping = methods::new("GeneGOMapping"), mapping = new.env(parent=emptyenv())),
					  validity = function(object) validObject(object@ggmapping)
		 )

setMethod("initialize", "GSminer", function(.Object, inputFile, sep, ...)
		  {
		  	.Object <- methods::callNextMethod(.Object, ...)
		  	.Object@ggmapping<-methods::new("GeneGOMapping", inputFile, sep, ...)
			mapping <- new.env(parent=emptyenv())
			mapping$BPoffsprings <- list()
			mapping$CCoffsprings <- list()
		  	.Object
		  }
)
