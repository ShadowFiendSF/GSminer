#include <R.h>
#include <Rinternals.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
extern SEXP sortpasteC(SEXP data);
SEXP sortpasteC(SEXP data)
{
    SEXP pasted;
    size_t i = 0;
    size_t n = length(VECTOR_ELT(data ,0));
    PROTECT(pasted = allocVector(STRSXP, n));
    SEXP col1 = VECTOR_ELT(data, 0);
    SEXP col2 = VECTOR_ELT(data, 1);
    for(; i < n; i++)
    {
        const char* str1 = CHAR(asChar(STRING_ELT(col1, i)));
        const char* str2 = CHAR(asChar(STRING_ELT(col2, i)));
        size_t len1 = strlen(str1);
        size_t len2 = strlen(str2);
        char* tmp = (char*)malloc(len1 + len2 + 1);
        memset(tmp, 0, len1 + len2 + 1);
        if(strcmp(str1, str2) >= 0)
        {
            strcat(tmp, str1);
            strcat(tmp, str2);
        }else{
            strcat(tmp, str2);
            strcat(tmp, str1);
        }
        SET_STRING_ELT(pasted, i, mkChar((const char*)tmp));
       free(tmp);
    }
    UNPROTECT(1);
    return pasted;
}
