// This file is prepared to register .call entry points.
// which was needed at the packaging step.
// For more information ->  "Writing R Extentions,  Registering native routines".

#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/RConverters.h>
#include <R_ext/Rdynload.h>

SEXP conductance_computation(SEXP,SEXP,SEXP);
SEXP  maximum_of_rows(SEXP);

static const R_CallMethodDef R_CallDef[] = {
	{"conductance_computation", (DL_FUNC)&conductance_computation, 4},
    {"maximum_of_rows", (DL_FUNC)&maximum_of_rows, 1},
    {NULL, NULL, 0},
};

void R_init_SamSPECTRAL(DllInfo *info)
{
  R_registerRoutines(info,NULL,R_CallDef,NULL,NULL);
}

