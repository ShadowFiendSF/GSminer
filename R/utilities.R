#' Format conversion
#'
#' @title convert to graph format
#' @param mapping the mapping which convert to igraph object
#' @return igraph object
#' @importFrom  AnnotationDbi toTable
#' @importFrom  GO.db GOBPPARENTS
#' @author Li Zhaohong && Wu Zefeng
#' @export
#' @examples
#' library(GO.db)
#' GO2igraph(mapping = GOBPPARENTS)
#'

GO2igraph <- function(mapping = GOBPPARENTS)
{
  if(! methods::is(mapping, "AnnDbBimap"))
    stop("the map object must be AnnDbBimap")
  BPtable <- toTable(mapping)
  graphBP <- igraph::graph.data.frame(BPtable)
}


#' Compute the GO term similarity
#' @title getGOSim
#' @param g igraph object
#' @param selectedGO list which name is selected GO term
#' @param method method for computing the similarity
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. 
#' @return similarity matrix 
#' @importFrom dnet dDAGreverse
#' @importFrom dnet dDAGannotate
#' @importFrom igraph V
#' @importFrom dnet dDAGtermSim
#' @export
#' @author Li Zhaohong && Wu Zefeng
#' @examples
#' 
#' library(igraph)
#' #just for an example:
#' #HPPA <-dRDataLoader(RData='ig.HPPA')
#' #g <- ig.HPPA
#' #getGOSim(g, selectedGO, method= "Resnik", multicores = 3)
#' 



getGOSim <- function(g, selectedGO, method=c("Resnik", "Lin", "Schlicker", "Jiang", "Pesquita"), multicores =NULL)
{
  g<-dDAGreverse(g)
  method <- match.arg(method)
  dag <- dDAGannotate(g, selectedGO, path.mode="all_paths",verbose=TRUE) 
 # data <- sapply(V(dag)$annotations, length)  
 # names(data) <- V(dag)$name
  sim <- dDAGtermSim(dag, terms = names(selectedGO), verbose = TRUE, method = method, multicores = multicores)
  sim <- as.matrix(sim)
  sim[!is.finite(sim)] <- 0
  sim
}

#' Clustering by apcluster to get the optimal cluster number and get the final go terms
#' @title GOCluster
#' @param sim similarity matrix 
#' @param selectedGO list which name is selected GO term
#' @param  seed for reproducibility,
#' @return final Selected GO terms
#' @importFrom apcluster apcluster
#' @importFrom dnet dDAGannotate
#' @importFrom igraph V
#' @importFrom dnet dDAGtermSim
#' @author Li Zhaohong && Wu Zefeng
#' @examples
#' \dontrun{
#' library(apcluster)
#' GOCluster(sim, selectedGO, seed = 100)
#' }
#'

GOCluster <- function(sim, selectedGO, seed)
{
  if(!is.matrix(sim)) sim <- as.matrix(sim)
  if(any(is.na(sim)) || any(is.infinite(sim)))
    stop("similarity matrix cannot contain infinite or NA.")
  if(dim(sim)[1] <= 1 || dim(sim)[2] <= 1)
    stop("similarity matrix must have at least one dimension.")
  finalSelect <- vector(mode = "character")
  apres <- apcluster(sim,seed=1000)
  message("GO cluster Number: ", length(apres))
  clusters <- apres@clusters
  names(clusters) <- names(apres@exemplars)
  meanSize <-round(mean(elemLen(selectedGO)))
  for (name in names(clusters)){
    if (length(selectedGO[name][[1]]) > meanSize){
      lenList <- c()
      for (n in names(clusters[name][[1]])){
        lenList <- append(lenList,length(selectedGO[n][[1]]))
      }
      elemCluster <- names(clusters[name][[1]])
      match <- lenList[which(abs(lenList - meanSize) == min(abs(lenList - meanSize)))][1]
      selectGO <- elemCluster[lenList == match][1] 
      finalSelect <- append(finalSelect, selectGO)
    }
    else {
      finalSelect <- append(finalSelect, name)
    }
  }
  finalSelect
}

#' @useDynLib GSminer elemLenC
elemLen <- function(l)
{
  if(!is.list(l)) return(NULL)
  .Call("elemLenC", l, PACKAGE="GSminer")
}

#' @useDynLib GSminer sortpasteC
sortpaste <- function(x)
{
	.Call("sortpasteC", x)
}







