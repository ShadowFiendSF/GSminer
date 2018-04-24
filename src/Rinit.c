#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

extern SEXP elemLenC(SEXP);
extern SEXP sortpasteC(SEXP);
//SEXP elemLen1(SEXP);

#if _MSC_VER >= 1000
__declspec(dllexport)
#endif


static const R_CallMethodDef R_CallDef[] = {
    {"elemLenC", (DL_FUNC) &elemLenC, 1},
    {"sortpasteC", (DL_FUNC) &sortpasteC, 1},
    {NULL, NULL, 0},
};

void R_init_GSminer(DllInfo *info)
{
  R_registerRoutines(info,NULL,R_CallDef,NULL,NULL);
  R_useDynamicSymbols(info, TRUE);
}
