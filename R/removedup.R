#' Remove duolicated rows
#' @param data data frame
#' @return an unique data frame
#for remove duolicated rows
#' @export
#' @examples
#' # GPSimNeg a data frame
#' test <- data.frame("a" = c("ab", "cd", "ef"), "b" = c("cd", "ab", "ef")) 
#' test <- removeDup(test)
#'


removeDup <- function(data)
{
	if(!identical(class(data), "data.frame"))
	{
		stop("in removeDup: input must be data frame!")
	}
	if(dim(data)[2] < 2)
	{
		stop("in removeDup: input must at least have tow columns!")
	}
	if(!is.character(data[,1])) data[,1] <- as.character(data[,1]) 
	if(!is.character(data[,2])) data[,2] <- as.character(data[,2]) 
	sp <- sortpaste(data)
	data[!duplicated(sp),]
}

