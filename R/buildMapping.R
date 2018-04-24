#' @include AllClass.R AllGenerics.R
NULL
#' @include AllClass.R AllGenerics.R
#' build Gene and GO mapping
#' @name buildMapping
#' @docType methods
#' @rdname buildMapping-methods
#' @title buildMapping method 
#' @param object the GeneGOMapping object
#' @param header input mapping file has a header or not
#' @return a GeneGOMapping with list of GO to Gene mapping
#' @exportMethod buildMapping
#' @aliases buildMapping,GeneGOMapping-method
#' @author Li Zhaohong && Wu Zefeng
#' @examples
#' #relationship between the gene and GO terms
#' #sep used in the input file by default is tab
#' inputFile <- paste0(system.file(package = "GSminer"), "/extdata/TAIR.GO")
#' miner1 <- GSminer(inputFile = inputFile, sep = "\t")
#' buildMapping(miner1@ggmapping)

setMethod("buildMapping", signature(object = "GeneGOMapping"),
		  function(object, header = FALSE)
		  {
				methods::validObject(object)
				gene2go <- utils::read.table(object@inputFile, sep = object@sep, header = header, stringsAsFactors = FALSE)
		  	  	gene2go <- gene2go[1:2]
		      	names(gene2go) <- c("gene", "go")
			  	go2gene <- tapply(gene2go$gene, as.factor(gene2go$go), 
			  														function(x) 
				  														as.vector(unique(x)), 
				  				simplify =FALSE)
#			    message("GO terms number: ", length(go2gene), " Gene number: ", length(unique(gene2go$go)))
				name <- dimnames(go2gene)[[1]]
		    	attr(go2gene, "dim") <- NULL
		    	attr(go2gene, "dimnames") <- NULL
		    	names(go2gene) <- name
			    assign("mapping", as.list(go2gene), envir = object@mapping)
		  }
)



#' gene for GO Biological Process (BP) terms 
#' @include AllClass.R AllGenerics.R
#' @name gene2BP
#' @docType methods
#' @rdname gene2BP-methods
#' @title gene2BP method 
#' @param object the GeneGOMapping object
#' @return a list of GO to Gene mapping
#' @exportMethod gene2BP
#' @importFrom GO.db GOBPOFFSPRING
#' @aliases gene2BP,GeneGOMapping-method
#' @author Li Zhaohong && Wu Zefeng
#' @examples
#'
#' inputFile <- paste0(system.file(package = "GSminer"), "/extdata/TAIR.GO")
#' #inputFile the name of input file which contain the mapping 
#' #relationship between the gene and GO terms
#' #sep used in the input file by default is tab
#' #for an example how to use:
#' miner <- GSminer(inputFile = inputFile, sep = "\t")
#' #gene2BP(miner@ggmapping)

setMethod("gene2BP", signature(object = "GeneGOMapping"),
		  function(object)
		  {
			methods::validObject(object)
			BP<-as.list(GOBPOFFSPRING)
			GO2gene <- vector("list",length(names(BP)))
            mapping <- NULL
            if (!is.null(r <- get0("mapping", envir = object@mapping, inherits = FALSE))) 
            {
                mapping <- r                
            }
		    if(is.null(mapping))
		    {
				buildMapping(object)
				mapping <- get("mapping", envir = object@mapping, inherits = FALSE)
		    }			
		  	names(GO2gene) <- names(BP)
		  	for (BP in names(GO2gene))
		  	{
				if(length(mapping[BP][[1]])!=0)
				{
				    GO2gene[BP][[1]] <- mapping[BP][[1]]
		    	}
			}
		  	message("  Biological prcess GO terms number: ", length(GO2gene), "   Gene number:  ", length(unique(unlist(GO2gene))))
			return(GO2gene)
		  }
)
#' @include AllClass.R AllGenerics.R
#' gene for Cytoplasm and nucleus Cell compartments(CC) GO terms
#' @name CC2gene
#' @docType methods
#' @rdname CC2gene-methods
#' @title CC2gene method 
#' @param object the GeneGOMapping object
#' @return a list of GO to Gene mapping
#' @exportMethod CC2gene
#' @importFrom GO.db GOBPOFFSPRING
#' @aliases CC2gene,GeneGOMapping-method
#' @author Li Zhaohong && Wu Zefeng
#' @examples
#'
#' inputFile <- paste0(system.file(package = "GSminer"), "/extdata/TAIR.GO")
#' miner <- GSminer(inputFile = inputFile, sep = "\t")
#' #for an example how to use:
#' #CC2gene(miner@ggmapping)
#'

setMethod("CC2gene", signature(object = "GeneGOMapping"),
		   function(object)
		   {
				CC <- as.list(GOCCOFFSPRING)
				CC2gene <- vector("list",length(names(CC)))
                mapping <- NULL
                if (!is.null(r <- get0("mapping", envir = object@mapping, inherits = FALSE))) 
                {
                    mapping <- r                
                }
				if(is.null(mapping))
				{
					buildMapping(object)
					mapping <- get("mapping", envir = object@mapping, inherits = FALSE)
				}
				names(CC2gene)<-names(CC)
                for (CC in names(CC2gene))
                {
					if(length(mapping[CC][[1]])!=0)
                    {
                		CC2gene[CC][[1]] <- mapping[CC][[1]]
                    }
		   		}
				message("Biological prcess GO terms number:", length(CC2gene), "Gene number: ", unique(unlist(CC2gene)))
				return(CC2gene)
		   }
)
