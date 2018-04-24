#' @include AllClass.R
NULL

#' Method buildMapping.
#
#' @name buildMapping
#' @rdname buildMapping-methods
# title buildMapping method 
# @exportMethod buildMapping
# author Li Zhaohong && Wu Zefeng
if ( !isGeneric("buildMapping") )
	setGeneric("buildMapping", function(object, ...) standardGeneric("buildMapping"))


#' Method gene2BP.
# docType methods
#' @name gene2BP
#' @rdname gene2BP-methods
# @exportMethod gene2BP
# title gene2BP method  
# author Li Zhaohong && Wu Zefeng
if ( !isGeneric("gene2BP") )
	setGeneric("gene2BP", function(object, ...) standardGeneric("gene2BP"))


#' Method CC2gene
# docType methods
#' @name CC2gene
#' @rdname CC2gene-methods
# title CC2gene method
# @exportMethod CC2gene
# author Li Zhaohong && Wu Zefeng
if ( !isGeneric("CC2gene") )
	setGeneric("CC2gene", function(object, ...) standardGeneric("CC2gene"))


#' loadCCoffsprings method
# docType methods
#' @name loadCCoffsprings
#' @rdname loadCCoffsprings-methods
# title loadCCoffsprings method  
# @exportMethod loadCCoffsprings
# author Li Zhaohong && Wu Zefeng
if( !isGeneric("loadCCoffsprings") )
	setGeneric("loadCCoffsprings", function(object, ...) standardGeneric("loadCCoffsprings"))


#' loadBPoffsprings method generics
# docType methods
#' @name loadBPoffsprings
#' @rdname loadBPoffsprings-methods
# title loadBPoffsprings method 
# @exportMethod loadBPoffsprings
# author Li Zhaohong && Wu Zefeng
if( !isGeneric("loadBPoffsprings") )
	setGeneric("loadBPoffsprings", function(object, ...) standardGeneric("loadBPoffsprings"))


#' getThreshold method generics
# docType methods
#' @name getThreshold
#' @rdname getThreshold-methods
# title getThreshold method
# @exportMethod getThreshold
# author Li Zhaohong && Wu Zefeng
if( !isGeneric("getThreshold") )
	setGeneric("getThreshold", function(object, ...) standardGeneric("getThreshold"))


#' selectGO method generics
# docType methods
#' @rdname selectGO-methods
#' @name selectGO
# title selectGO method 
# exportMethod selectGO
# author Li Zhaohong && Wu Zefeng
if( !isGeneric("selectGO") )
	setGeneric("selectGO", function(object, ...) standardGeneric("selectGO"))


#' outputPosNeg method generics
# docType methods
#' @rdname outputPosNeg-methods
#' @name outputPosNeg
# title outputPosNeg method
# @exportMethod outputPosNeg
# author Li Zhaohong && Wu Zefeng
if( !isGeneric("outputPosNeg") )
	setGeneric("outputPosNeg", function(object, ...) standardGeneric("outputPosNeg"))

#' filterNeg method generics
# docType methods
#' @name filterNeg
#' @rdname filterNeg-methods
# @exportMethod filterNeg
# title filterNeg method  
# author Li Zhaohong && Wu Zefeng
if( !isGeneric("filterNeg") )
	setGeneric("filterNeg", function(object, ...) standardGeneric("filterNeg"))
