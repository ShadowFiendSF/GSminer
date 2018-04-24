#include <R.h>
#include <Rinternals.h>
#include <ctype.h>

SEXP elemLenC(SEXP x)
{
    size_t n = length(x);
    if(!Rf_isNewList(x))
        error("int elemlen expected type %s\n", type2char(TYPEOF(x)));
    SEXP ans;
    PROTECT(ans = allocVector(INTSXP, n));
    memset(INTEGER(ans), 0, n * sizeof(int));
    for(int i = 0; i < n; i++)
    {
	if(Rf_isNull(VECTOR_ELT(x, i)))
        {
	   INTEGER(ans)[i] = 0;
	}else{
	   INTEGER(ans)[i] = length(VECTOR_ELT(x, i));
	} 
    }
    setAttrib(ans, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    UNPROTECT(1);
    return(ans);
}

