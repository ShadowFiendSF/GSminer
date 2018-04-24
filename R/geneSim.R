# 
#
#' The similarity between genes(modified the dDAGgeneSim in dnet)
#' @importFrom dnet dDAGinduce
#' @importFrom dnet dDAGtip
#' @importFrom dnet dCheckParallel
#' @importFrom foreach foreach
#' @import igraph
#' @param g  an object of class "igraph" 
#' @param genes the genes between which pair-wise semantic similarity is calculated
#' @param method.gene  the method used for how to derive semantic similarity
#'                     between genes from semantic similarity between terms
#' @param method.term the method used to measure semantic similarity between terms
#' @param force  logical 
#' @param fast logical to indicate whether a vectorised fast computation is used
#' @param parallel logical to indicate whether parallel computation
#' @param verbose logical
#' @param multicores parallel computation
#' @return  a sparse matrix of similarity between input genes
#' @export
#' @examples
#' #load HPPA as igraph object
#' HPPA <-dRDataLoader(RData='ig.HPPA')
#' g <- ig.HPPA
#' #load human genes annotated by HPPA
#' Hs.egHPPA <- dRDataLoader(RData='org.Hs.egHPPA')
#' dag <- dDAGannotate(g, annotations=org.Hs.egHPPA,
#' path.mode="all_paths", verbose=TRUE)
#' allgenes <- unique(unlist(V(dag)$annotations))
#' genes <- sample(allgenes,5)
#' #just for an example:
#' #geneSim(g = dag, genes = genes, method.gene="BM.average", 
#' #        method.term = "Resnik", parallel = TRUE, 
#' #        multicores = 2, verbose = FALSE)
#'
#'




geneSim <- function (g, genes=NULL, method.gene=c("BM.average","BM.max","BM.complete","average","max"), method.term=c("Resnik","Lin","Schlicker","Jiang","Pesquita"), 
	force=TRUE, fast=TRUE, parallel=TRUE, multicores=NULL, verbose=TRUE)
{  
    startT <- Sys.time()
    if(verbose){
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    message("", appendLF=TRUE)
    }
  ####################################################################################
  
    method.gene <- match.arg(method.gene)
    method.term <- match.arg(method.term)
  
    ig <- g
    if (class(ig) != "igraph"){
    stop("The function must apply to either 'igraph' or 'graphNEL' object.\n")
    }
  
    if(is.null(V(ig)$annotations) | is.null(V(ig)$IC)){
		stop("The function requires that input graph has already contained annotation data. Please first run 'dDAGannotate'.\n")
    }
  
  ####################################################
  ## A function to indicate the running progress
	progress_indicate <- function(i, B, step, flag=FALSE){
		if(i %% ceiling(B/step) == 0 | i==B | i==1){
      		if(flag & verbose){
       			 message(sprintf("\t%d out of %d (%s)", i, B, as.character(Sys.time())), appendLF=TRUE)
     		 }
    	}
	}
  
  ####################################################
  
  if(verbose){
    message(sprintf("First, extract all annotatable genes (%s)...", as.character(Sys.time())), appendLF=TRUE)
  }
  
  anno <- V(ig)$annotations
  allgenes <- sort(as.character(unique(unlist(anno))))
  
  ## checking input genes
  genes <- genes[!is.na(genes)]
  if(is.null(genes) || is.na(genes)){
    genes <- allgenes
  }else{
    flag <- genes %in% allgenes
    if(sum(flag)!=0){
      genes <- genes[flag]
    }else{
      genes <- allgenes
    }
  }
  
  ## pre-compute a sparse matrix of input genes x terms
  allterms <- 1:length(anno)
  sGT <- Matrix::Matrix(0, nrow=length(genes), ncol=length(allterms), sparse=TRUE)
  for(j in 1:length(allterms)){
    ind <- match(anno[[j]], genes)
    flag <- ind[!is.na(ind)]
    if(length(flag)!=0){
      sGT[flag,j] <- 1
    }
  }
  colnames(sGT) <- V(ig)$name
  rownames(sGT) <- genes
  
  if(verbose){
    message(sprintf("\tthere are %d input genes amongst %d annotatable genes", length(genes), length(allgenes)), appendLF=TRUE)
  }
  
  ## a list of genes, each containing terms annotated by
  genes2terms <- sapply(1:length(genes), function(x){
    res <- names(which(sGT[x,]==1))
    if(force){
      subg <- dDAGinduce(ig, nodes_query=res, path.mode="all_paths")
      res <- dDAGtip(subg)
    }
    return(res)
  })
  names(genes2terms) <- genes
  terms <- unique(unlist(genes2terms))
  
  ## also instore index for terms (in genes2terms)
  genes2terms_index <- sapply(genes2terms, function(x){
    match(x, terms)
  })
  
  if(verbose){
    if(force){
      message(sprintf("Second, pre-compute semantic similarity between %d terms (forced to be the most specific for each gene) using %s method (%s)...", length(terms), method.term, as.character(Sys.time())), appendLF=TRUE)
    }else{
      message(sprintf("Second, pre-compute semantic similarity between %d terms using %s method (%s)...", length(terms), method.term, as.character(Sys.time())), appendLF=TRUE)
    }
  }
  ## pre-compute semantic similarity between terms in subject
  sim.term <- suppressMessages(dDAGtermSim(ig, terms=terms, method=method.term, parallel=parallel, multicores=multicores, verbose=TRUE))
  
  if(verbose){
    message(sprintf("Last, calculate pair-wise semantic similarity between %d genes using %s method (%s)...", length(genes), method.gene, as.character(Sys.time())), appendLF=TRUE)
  }
  num_genes <- length(genes2terms)
  
  ###### parallel computing
  flag_parallel <- FALSE
  if(parallel==TRUE){
    flag_parallel <- dCheckParallel(multicores=multicores, verbose=verbose)
    if(flag_parallel){
      if(method.gene=='average'){
        i <- 1
        sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_genes-1), .inorder=TRUE, .combine=rbind), {
          ind1 <- genes2terms_index[[i]]
          progress_indicate(i, num_genes, 10, flag=TRUE)
          fast <- TRUE
          if(fast){
            js <- (i+1):num_genes
            ind_js <- genes2terms_index[js]
            sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
            new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
            res <- sapply(1:length(ind_js), function(k){
              mean(sim12[,which(new_ind_js==k)])
            })
            x <- rep(0, num_genes)
            x[js] <- res
            x
          }
        })
      }else if(method.gene=='max'){
        i <- 1
        sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_genes-1), .inorder=TRUE, .combine=rbind), {
          ind1 <- genes2terms_index[[i]]
          progress_indicate(i, num_genes, 10, flag=TRUE)
          fast <- TRUE
          if(fast){
            js <- (i+1):num_genes
            ind_js <- genes2terms_index[js]
            sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
            new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
            res <- sapply(1:length(ind_js), function(k){
              max(sim12[,which(new_ind_js==k)])
            })
            x <- rep(0, num_genes)
            x[js] <- res
            x
          }
        })
      }else if(method.gene=='BM.average'){
        i <- 1
        sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_genes-1), .inorder=TRUE, .combine=rbind), {
          ind1 <- genes2terms_index[[i]]
          progress_indicate(i, num_genes, 10, flag=TRUE)
          fast <- TRUE
          if(fast){
            js <- (i+1):num_genes
            ind_js <- genes2terms_index[js]
            sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
            new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
            res <- sapply(1:length(ind_js), function(k){
              x <- as.matrix(sim12[,which(new_ind_js==k)])
              0.5*(mean(apply(x,1,max)) + mean(apply(x,2,max)))
            })
            x <- rep(0, num_genes)
            x[js] <- res
            x
          }
        })
      }else if(method.gene=='BM.max'){
        i <- 1
        sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_genes-1), .inorder=TRUE, .combine=rbind), {
          ind1 <- genes2terms_index[[i]]
          progress_indicate(i, num_genes, 10, flag=TRUE)
          fast <- TRUE
          if(fast){
            js <- (i+1):num_genes
            ind_js <- genes2terms_index[js]
            sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
            new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
            res <- sapply(1:length(ind_js), function(k){
              x <- as.matrix(sim12[,which(new_ind_js==k)])
              max(mean(apply(x,1,max)), mean(apply(x,2,max)))
            })
            x <- rep(0, num_genes)
            x[js] <- res
            x
          }
        })
      }else if(method.gene=='BM.complete'){
        i <- 1
        sim <- foreach::`%dopar%` (foreach::foreach(i=1:(num_genes-1), .inorder=TRUE, .combine=rbind), {
          ind1 <- genes2terms_index[[i]]
          progress_indicate(i, num_genes, 10, flag=TRUE)
          fast <- TRUE
          if(fast){
            js <- (i+1):num_genes
            ind_js <- genes2terms_index[js]
            sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
            new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
            res <- sapply(1:length(ind_js), function(k){
              x <- as.matrix(sim12[,which(new_ind_js==k)])
              min(c(apply(x,1,max),apply(x,2,max)))
            })
            x <- rep(0, num_genes)
            x[js] <- res
            x
          }
        })
      }
      
      ## add the last row
      sim <- rbind(sim, rep(0, num_genes))
      
      sim <- sim + Matrix::t(sim)
      sim <- Matrix::Matrix(sim, sparse=TRUE)
    }
  }
  
  ###### non-parallel computing
  if(flag_parallel==FALSE){
    ## calculate pair-wise semantic similarity between input genes
    sim <- Matrix::Matrix(0, nrow=length(genes), ncol=length(genes), sparse=TRUE)
    
    ## print with possibly greater accuracy:
    ##op <- options(digits.secs = 6)
    ##options(op)
    
    if(method.gene=='average'){
      for(i in 1:(num_genes-1)){
        ind1 <- genes2terms_index[[i]]
        progress_indicate(i, num_genes, 10, flag=TRUE)
        if(fast){
          js <- (i+1):num_genes
          ind_js <- genes2terms_index[js]
          sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
          new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
          res <- sapply(1:length(ind_js), function(k){
            mean(sim12[,which(new_ind_js==k)])
          })
          sim[i,js] <- res
        }else{
          for(j in (i+1):num_genes){
            ind2 <- genes2terms_index[[j]]
            ## pairwise similarity between terms
            sim12 <- as.matrix(sim.term[ind1, ind2])
            sim[i,j] <- mean(sim12)
          }
        }
      }
    }else if(method.gene=='max'){
      for(i in 1:(num_genes-1)){
        ind1 <- genes2terms_index[[i]]
        progress_indicate(i, num_genes, 10, flag=TRUE)
        if(fast){
          js <- (i+1):num_genes
          ind_js <- genes2terms_index[js]
          sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
          new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
          res <- sapply(1:length(ind_js), function(k){
            max(sim12[,which(new_ind_js==k)])
          })
          sim[i,js] <- res
        }else{
          for(j in (i+1):num_genes){
            ind2 <- genes2terms_index[[j]]
            ## pairwise similarity between terms
            sim12 <- as.matrix(sim.term[ind1, ind2])
            sim[i,j] <- max(sim12)
          }
        }
      }
    }else if(method.gene=='BM.average'){
      for(i in 1:(num_genes-1)){
        ind1 <- genes2terms_index[[i]]
        progress_indicate(i, num_genes, 10, flag=TRUE)
        if(fast){
          js <- (i+1):num_genes
          ind_js <- genes2terms_index[js]
          sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
          new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
          res <- sapply(1:length(ind_js), function(k){
            x <- as.matrix(sim12[,which(new_ind_js==k)])
            0.5*(mean(apply(x,1,max)) + mean(apply(x,2,max)))
          })
          sim[i,js] <- res
        }else{
          for(j in (i+1):num_genes){
            ind2 <- genes2terms_index[[j]]
            ## pairwise similarity between terms
            sim12 <- as.matrix(sim.term[ind1, ind2])
            sim[i,j] <- 0.5*(mean(apply(sim12,1,max)) + mean(apply(sim12,2,max)))
          }
        }
      }
      
    }else if(method.gene=='BM.max'){
      for(i in 1:(num_genes-1)){
        ind1 <- genes2terms_index[[i]]
        progress_indicate(i, num_genes, 10, flag=TRUE)
        if(fast){
          js <- (i+1):num_genes
          ind_js <- genes2terms_index[js]
          sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
          new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
          res <- sapply(1:length(ind_js), function(k){
            x <- as.matrix(sim12[,which(new_ind_js==k)])
            max(mean(apply(x,1,max)), mean(apply(x,2,max)))
          })
          sim[i,js] <- res
        }else{
          for(j in (i+1):num_genes){
            ind2 <- genes2terms_index[[j]]
            ## pairwise similarity between terms
            sim12 <- as.matrix(sim.term[ind1, ind2])
            sim[i,j] <- max(mean(apply(sim12,1,max)), mean(apply(sim12,2,max)))
          }
        }
      }
    }else if(method.gene=='BM.complete'){
      for(i in 1:(num_genes-1)){
        ind1 <- genes2terms_index[[i]]
        progress_indicate(i, num_genes, 10, flag=TRUE)
        if(fast){
          js <- (i+1):num_genes
          ind_js <- genes2terms_index[js]
          sim12 <- matrix(sim.term[ind1, unlist(ind_js)], nrow=length(ind1))
          new_ind_js <- rep(1:length(ind_js), sapply(ind_js,length))
          res <- sapply(1:length(ind_js), function(k){
            x <- as.matrix(sim12[,which(new_ind_js==k)])
            min(c(apply(x,1,max),apply(x,2,max)))
          })
          sim[i,js] <- res
        }else{
          for(j in (i+1):num_genes){
            ind2 <- genes2terms_index[[j]]
            ## pairwise similarity between terms
            sim12 <- as.matrix(sim.term[ind1, ind2])
            sim[i,j] <- min(c(apply(sim12,1,max),apply(sim12,2,max)))
          }
        }
      }
    }
    sim <- sim + Matrix::t(sim)
    
  }
  
  rownames(sim) <- colnames(sim) <- genes
  
  ####################################################################################
  endT <- Sys.time()
  if(verbose){
    message("", appendLF=TRUE)
    message(paste(c("Finish at ",as.character(endT)), collapse=""), appendLF=TRUE)
  }
  
  runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
  message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
  
  invisible(sim)
}
