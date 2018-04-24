#' @include AllClass.R AllGenerics.R
NULL

#' Output negative gene pairs
#' @name outputPosNeg
#' @docType methods
#' @rdname outputPosNeg-methods
#' @title outputPosNeg method 
#' @param u the upper bound, 1/2 by default 
#' @param l the lower bound, 1/4 by default
#' @param object the GSminer object
#' @param method.term method for computing the similarity
#' @param method.gene the method used for how to derive semantic similarity between genes from semantic similarity between terms.
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. 
#' @param seed for reproducibility, the seed of the random number generator can be set to a fixed value before adding noise.
#' @param posFilename postive gene output file name
#' @param negFilename negative gene output file name without filtered by CC GO terms
#' @param filterByCC weather filtered by CC GO terms by default is TRUE
#' @param FNegFile negative gene output file name filtered by CC GO terms
#' @param verbose logical to indicate whether the messages will be displayed in the screen.
#' @param sep the delimiter of the file
#' @return a vector which name is selected GO term
#' @importFrom AnnotationDbi toTable
#' @importFrom GO.db GOBPPARENTS
#' @importFrom dplyr group_by
#' @importFrom dplyr filter
#' @importFrom dplyr ungroup
#' @importFrom dplyr select
#' @importFrom dplyr n
#' @importFrom rlang sym
#' @import magrittr
#' @import utils
#' @exportMethod outputPosNeg
#' @aliases outputPosNeg,GSminer-method
#' @author Lizhaohong && Wu Zefeng
#' @examples
#'
#' inputFile <- paste0(system.file(package = "GSminer"), "/extdata/TAIR.GO")
#' miner <- GSminer(inputFile = inputFile, sep = "\t")
#' #outputPosNeg(miner, u = 0.1, l = 0.2, method.term = "Resnik", 
#' #             method.gene = "BM.average", multicores = 2, 
#' #             verbose = FALSE, posFilename = "pos.txt", 
#' #            filterByCC = FALSE, FNegFile = "neg.text")
#'


setMethod("outputPosNeg", signature(object = "GSminer"),
			function(object, u = 0.5, l = 0.25, method.term = c("Resnik", "Lin", "Schlicker", "Jiang", "Pesquita"), method.gene = c("BM.average", "BM.max", "BM.complete", "average", "max"),
                     multicores = NULL, verbose = TRUE, seed = NA, posFilename = "positiveGene.txt", negFilename = "negativeGene.txt", sep = "\t", filterByCC = TRUE, FNegFile = "negativeGene.filtered.txt")
			{
                if(!exists("BPoffsprings", envir = object@mapping))
                {
                    loadBPoffsprings(object)
                }
                options(stringsAsFactors=FALSE)
                method.term <- match.arg(method.term)
                method.gene <- match.arg(method.gene)
				BPoffsprings <- get("BPoffsprings", envir = object@mapping)
                finalSelectedGO <- selectGO(object, u = u, l = l, method = method.term, multicores = multicores, seed = seed)
                message("final Selected GO Number: ", length(finalSelectedGO))
                finalSelectedGene <- BPoffsprings[finalSelectedGO]
                annotations <- BPoffsprings[finalSelectedGO]
                dag <- dDAGannotate(dDAGreverse(GO2igraph(mapping = GOBPPARENTS)), annotations, path.mode="all_paths",verbose=TRUE)
                message("\n\ncompute the similarity between genes:\n")
                allAnnotatedGenesCorr <- geneSim(g = dag, genes = unique(unlist(finalSelectedGene)), 
                                         method.gene = method.gene, method.term = method.term, parallel = TRUE, 
                                         multicores = multicores, verbose = verbose)
                allAnnotatedGenesCorr <- as.matrix(allAnnotatedGenesCorr)
                GPSimPos <- data.frame()
                message("Generating the positive gene pairs......")
                #for (name in names(finalSelectedGene))
                #{
                #    temp <- finalSelectedGene[[name]]
                #    subMatrix <- allAnnotatedGenesCorr[temp,temp]
                #    GPSimSub <- data.frame(row = rownames(subMatrix)[row(subMatrix)[upper.tri(subMatrix)]],
                #    col = colnames(subMatrix)[col(subMatrix)[upper.tri(subMatrix)]], corr = subMatrix[upper.tri(subMatrix)], class = "pos")
                #   GPSimPos <- rbind(GPSimPos, GPSimSub)
                #}
                for (name in names(finalSelectedGene))
                {
                    temp <- finalSelectedGene[[name]]
                    if(length(temp) == 0)
                        next
                    tempdf <- expand.grid(temp, temp)
                    tempdf <- tempdf[!tempdf[,1] == tempdf[,2],]
                    GPSimPos <- rbind(GPSimPos, tempdf)
                }
                corrP <- vector("numeric", dim(GPSimPos)[1])
                indxP1 <- match(GPSimPos[,1], rownames(allAnnotatedGenesCorr))
                indxP2 <- match(GPSimPos[,2], colnames(allAnnotatedGenesCorr))
                for(i in 1:length(indxP1))
                {
                    corrP[i] <- allAnnotatedGenesCorr[indxP1[i], indxP2[i]]               
                }
                GPSimPos[, 3] <- corrP
                message("Remove the duplicated  positive gene pairs......") 
  
                #GPSimPos <- GPSimPos[!duplicated(t(apply(GPSimPos, 1, sort))),]
                names(GPSimPos) <- c("gene1", "gene2", "Sim")
                GPSimPos <- removeDup(GPSimPos)
                # CODE SNIPPET used for select the duplicated rows but it did not used here                 
                #GPSimPos <- GPSimPos %>%
                #            group_by(pmin(gene1, gene2), pmax(gene1, gene2)) %>%
                #            filter(n() >= 2) %>%
                #            ungroup() %>%
                #           select(gene1, gene2, Sim, Class)

                message("Output the positive gene pairs......")
                write.table(GPSimPos, file = posFilename, quote = FALSE, append = FALSE, row.names = FALSE, col.names = FALSE,sep= sep)
                message("Generating the negative gene pairs......")
                positiveGoTerms <- data.frame(t(combn(finalSelectedGO, m = 2)),stringsAsFactors = FALSE)

                num <- 0
                for (r in 1:nrow(positiveGoTerms))
                {
                    set1 <- finalSelectedGene[[positiveGoTerms[ r, 1 ]]]
                    set2 <- finalSelectedGene[[positiveGoTerms[ r, 2 ]]]
                    num <- num + length(set1) * length(set2)
                }

                gene1 <- vector("character", num)
                gene2 <- vector("character", num)
                corr <- vector("numeric", num)
                marker <- 1                
                for (r in 1:nrow(positiveGoTerms))
                {                   
                    set1 <- finalSelectedGene[[positiveGoTerms[ r, 1 ]]]
                    set2 <- finalSelectedGene[[positiveGoTerms[ r, 2 ]]]
                    #subMatrix <- allAnnotatedGenesCorr[set1, set2]
                    #GPsim  <-as.data.frame(as.table(subMatrix))
                    cb <- expand.grid(set1, set2,  stringsAsFactors = FALSE)
                    gene1[marker : (marker + dim(cb)[1] - 1)] <- cb[,1]
                    gene2[marker : (marker + dim(cb)[1] - 1)] <- cb[,2]
                    marker <- marker + dim(cb)[1] 
                }
                indx1 <- match(gene1, rownames(allAnnotatedGenesCorr))
                indx2 <- match(gene2, colnames(allAnnotatedGenesCorr))
                for(i in 1:length(indx1))
                {
                    corr[i] <- allAnnotatedGenesCorr[indx1[i], indx2[i]]               
                }
                GPSimNeg <- data.frame("gene1" = gene1, "gene2" = gene2, "Corr" = corr, stringsAsFactors = FALSE)
             
                                
                message("Remove the duplicated  negative gene pairs......")
                #names(GPSimPos) <- c("gene1", "gene2", "Sim")
                GPSimNeg <- removeDup(GPSimNeg)
                #used for select the duplicated rows but it did not used here                
                #GPSimNeg <- GPSimNeg %>%
                #            group_by(pmin(gene1, gene2), pmax(gene1, gene2)) %>%
                #           filter(n() >= 2) %>%
                #            ungroup() %>%
                #            select(gene1, gene2, Sim)
                message("Output the negative gene pairs......")
                write.table(GPSimNeg, file = negFilename, quote = FALSE, append = FALSE, row.names = FALSE, col.names = FALSE, sep= sep)
                if(filterByCC == TRUE && dim(GPSimNeg)[1] != 0)
                {
                    filterNeg(object, negFilename, FNegFile, sep)
                }
			}
)


#' Filter neagtive paris by cytoplasm and nucleus cell compartments
#' @param object the GSminer object
#' @param negFilename the negative gene file name used for filtering
#' @param FNegFile the output file name which is filtered by cytoplasm and nucleus cell compartments
#' @param sep the delimiter of the file
#' @return neagtive paris file
#' @importFrom GO.db GOCCOFFSPRING
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr filter
#' @importFrom dplyr ungroup
#' @importFrom dplyr select
#' @importFrom dplyr n
#' @importFrom rlang sym
#' @importFrom rlang ":="
#' @importFrom rlang "!!"
#' @importFrom magrittr "%>%"
#' @import utils
#' @aliases filterNeg,GSminer-method
#' @examples
#'
#' \dontrun{
#' inputFile <- paste0(system.file(package = "GSminer"), "/extdata/TAIR.GO")
#' library(GSminer)
#' test <- GSminer(inputFile = inputFile)
#' filterNeg(object, negFilename = "negFilename", FNegFile = "FNegFile")
#' \}


setMethod("filterNeg", signature(object = "GSminer"),
            function(object, negFilename = "negativeGene.txt", FNegFile = "negativeGene.filtered.txt", sep = "\t")
            {
                message("Filter neagtive paris by cytoplasm and nucleus cell compartments!")
                negPairs <- utils::read.table( negFilename, header = FALSE, sep = sep, stringsAsFactors = FALSE)
                if(!exists("CCoffsprings", envir = object@mapping))
                {
                    loadCCoffsprings(object)
                }
                CCoffsprings <- get("CCoffsprings", envir = object@mapping)
                geneChl <- setdiff(CCoffsprings["GO:0005737"][[1]], CCoffsprings["GO:0005634"][[1]]) 
                geneNuc <- setdiff(CCoffsprings["GO:0005634"][[1]], CCoffsprings["GO:0005737"][[1]])
                in1 <- sym(sprintf("in1"))
                in2 <- sym(sprintf("in2"))
                in3 <- sym(sprintf("in3"))
                in4 <- sym(sprintf("in4"))
                V1 <- sym(sprintf("V1"))
                V2 <- sym(sprintf("V2"))
                result <- negPairs %>% 
                                mutate( !!in1 := (!!V1) %in% geneChl,
                                        !!in2 := (!!V1) %in% geneNuc, 
                                        !!in3 := (!!V2) %in% geneChl, 
                                        !!in4 := (!!V2) %in% geneNuc) %>%  
                                filter(((!!in1) & (!!in4)) 
                                      |((!!in2) & (!!in3)))

                result <- removeDup(result)                
                #result <- result[!duplicated(t(apply(result, 1, sort))),]
                #names(result) <- c("gene1", "gene2", "Sim")
                #result <- result %>%
                #          group_by(pmin(gene1, gene2), pmax(gene1, gene2)) %>%
                #          filter(n() >= 2) %>%
                #          ungroup() %>%
                #          select(gene1, gene2, Sim)              
                write.table(result, file = FNegFile, col.names = FALSE, quote = FALSE, sep = "\t", 
                             row.names = FALSE,append = TRUE)
            }

)
