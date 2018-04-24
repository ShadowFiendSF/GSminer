#' @include AllClass.R AllGenerics.R
NULL

#' select the GO term size
#' @name getThreshold
#' @docType methods
#' @rdname getThreshold-methods
#' @title getThreshold method 
#' @param object the GSminer object
#' @param u the upper bound, 1/2 by default 
#' @param l the lower bound, 1/4 by default
#' @return a list of upper bound and lower bound
#' @exportMethod getThreshold
#' @aliases getThreshold,GSminer-method
#' @author Li Zhaohong && Wu Zefeng
#' @examples
#'
#' inputFile <- paste0(system.file(package = "GSminer"), "/extdata/TAIR.GO")
#' miner <- GSminer(inputFile = inputFile, sep = "\t")
#' # just for an example:
#' #getThreshold(miner, u = 0.5, l = 0.25)


setMethod("getThreshold", signature(object = "GSminer"),
		  function(object, u = 0.5, l = 0.25)
		  {
			loadBPoffsprings(object)
			BPoffsprings <- get("BPoffsprings", envir = object@mapping)
			total <- length(unique(unlist(BPoffsprings)))
		  	result <- vector("list", length=2)
			names(result) <- c("upper", "lower")
			upper <- round(total*u)
			lower <- round(total*l)
			totalGene <- function(GOList, size)
			{
				if(!is.list(GOList)) return(NULL)
				length(
					   unique(
							  unlist(
									 GOList[elemLen(GOList) <= size]
									 )
							  )
					   )
			}
			message("Now selectting the upper threshold, please wait!")
			for(x in seq_along( 1:total))
			{
#		    	message("Now x: ", x, " ", totalGene(BPoffsprings, x), " >= ", upper)
				if(totalGene(BPoffsprings, x) >= upper)
				{
					result[["upper"]] <- x
					break
				 }
			 }
			message("Now selectting the lower threshold, please wait!")
		  	for (x in seq_along( 1:total))
			{
#		    	message("Now x: ", x, " ", totalGene(BPoffsprings, x), " >= ", lower)
				if(totalGene(BPoffsprings, x) >= lower)
				{
					result[["lower"]] <- x
					break
				}
			}
			if(result[["upper"]] == result[["lower"]])
			{
				warning("upper threshold and lower threshold has the same value. Please reset!")
			}		  
		  	return(result)
		  }
)

#' select the GO term based on the threshold
#' @name selectGO
#' @docType methods
#' @rdname selectGO-methods
#' @title selectGO method 
#' @param u the upper bound, 1/2 by default 
#' @param l the lower bound, 1/4 by default
#' @param object the GSminer object
#' @param method method for computing the similarity
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. 
#' @param seed for reproducibility, the seed of the random number generator can be set to a fixed value before adding noise.
#' @return a vector which name is selected GO term
#' @exportMethod selectGO
#' @aliases selectGO,GSminer-method
#' @author Li Zhaohong && Wu Zefeng
#' @examples
#'
#' inputFile <- paste0(system.file(package = "GSminer"), "/extdata/TAIR.GO")
#' miner <- GSminer(inputFile = inputFile, sep = "\t")
#' #just for an example
#' #selectGO(miner,u = 0.5, l = 0.25, method= "Resnik", multicores = 2)
#' 

setMethod("selectGO", signature(object = "GSminer"),
			function(object, u = 0.5, l = 0.25, method=c("Resnik", "Lin", "Schlicker", "Jiang", "Pesquita"), multicores = NULL, seed = NA)
			{
				method <- match.arg(method)
				threshold <- getThreshold(object, u, l)
				`%ni%` <-Negate("%in%")
				if("upper" %ni% names(threshold) && "lower" %ni% names(threshold))
					stop("threshold must have the upper bound and the lower bound")
				BPoffsprings <- list()
				if(!exists("BPoffsprings", envir = object@mapping))
				{
					loadBPoffsprings(object)
					assign("BPoffsprings", get("BPoffsprings", envir = object@mapping))
				}else{
					assign("BPoffsprings", get("BPoffsprings", envir = object@mapping))
				}
				selectedGO <- BPoffsprings[ elemLen(BPoffsprings) >= threshold[["lower"]] & elemLen(BPoffsprings) <= threshold[["upper"]] ]
				simMat <- getGOSim( GO2igraph(mapping = GOBPPARENTS), 
									selectedGO = selectedGO, 
									method = method, 
									multicores = multicores)
				GOCluster(sim = simMat, selectedGO = selectedGO, seed = seed)
			}
)

