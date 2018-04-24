#' @include AllClass.R AllGenerics.R
NULL

#' load Biological Process GO offsprings for genes
#' @name loadBPoffsprings
#' @docType methods
#' @rdname loadBPoffsprings-methods
#' @title loadBPoffsprings method 
#' @param object the GeneGOMapping object
#' @return a list of BP GO to Gene mapping
#' @importFrom GO.db GOBPOFFSPRING
#' @exportMethod loadBPoffsprings
#' @aliases loadBPoffsprings,GSminer-method
#' @author Li Zhaohong && Wu Zefeng
#' @examples
#' inputFile <- paste0(system.file(package = "GSminer"), "/extdata/TAIR.GO")
#' test <- GSminer(inputFile = inputFile)
#' #for an example:
#' #loadBPoffsprings(test)


setMethod("loadBPoffsprings", signature(object = "GSminer"),
		   function(object)
		   {
				BP <- as.list(GOBPOFFSPRING)
				BPoffsprings <- vector("list",length(names(BP)))
				names(BPoffsprings) <- names(BP)
				mapping <- list()
				if(!exists("mapping", envir = object@ggmapping@mapping) || length(mapping) == 0)
				{
					buildMapping(object@ggmapping)
					assign("mapping", get("mapping", envir = object@ggmapping@mapping, inherits = FALSE))
				}else{
					assign("mapping", get("mapping", envir = object@ggmapping@mapping, inherits = FALSE))
				}
				for (n in names(BPoffsprings))
				{
			    	temp<-c()
			    	for (go in BP[n][[1]])
					{
				    	temp <- append(temp, mapping[go][[1]]) 
				   	}
					temp <- append(temp, mapping[n][[1]])
					if(length(temp)!=0)
					{
						BPoffsprings[n][[1]]<-unique(temp)
					}
			   	}
				assign("BPoffsprings", BPoffsprings, envir = object@mapping)
			}
)

#' load cytoplasm and nucleus cell compartments GO terms for genes
#' @name loadCCoffsprings
#' @docType methods
#' @rdname loadCCoffsprings-methods
#' @title loadCCoffsprings method 
#' @param object the GeneGOMapping object
#' @return a list of CC GO to Gene mapping
#' @importFrom GO.db GOCCOFFSPRING
#' @exportMethod loadCCoffsprings
#' @aliases loadCCoffsprings,GSminer-method
#' @author Li Zhaohong && Wu Zefeng
#' @examples
#' inputFile <- paste0(system.file(package = "GSminer"), "/extdata/TAIR.GO")
#' test <- GSminer(inputFile = inputFile)
#' #just like loadBPoffsprings
#' #loadCCoffsprings(test)

setMethod("loadCCoffsprings", signature(object = "GSminer"),
		   function(object)
		   {
		    	CC<-as.list(GOCCOFFSPRING)
				CCoffsprings <- vector("list",length(names(CC)))
				names(CCoffsprings) <- names(CC)
				mapping <- list()
				if(!exists("mapping", envir = object@ggmapping@mapping) || length(mapping) == 0)
				{
					buildMapping(object@ggmapping)
					assign("mapping", get("mapping", envir = object@ggmapping@mapping, inherits = FALSE))
				}else{
					assign("mapping", get("mapping", envir = object@ggmapping@mapping, inherits = FALSE))
				}
				for (n in names(CCoffsprings))
				{
					temp <- c()
					for (go in CC[n][[1]])
					{
						temp <- append(temp, mapping[go][[1]]) 
					}
					temp <- append(temp, mapping[n][[1]])
					if(length(temp)!=0)
					{
						CCoffsprings[n][[1]]<-unique(temp)
					}
				}
				assign("CCoffsprings", CCoffsprings, envir = object@mapping)
		   }
)
