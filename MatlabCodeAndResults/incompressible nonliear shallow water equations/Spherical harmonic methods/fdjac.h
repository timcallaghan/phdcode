/*
 * MATLAB Compiler: 2.1
 * Date: Tue Mar 09 12:17:15 2004
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C"
 * "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "fdjac" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __fdjac_h
#define __fdjac_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_fdjac(void);
extern void TerminateModule_fdjac(void);
extern _mexLocalFunctionTable _local_function_table_fdjac;

extern mxArray * mlfFdjac(mxArray * x, mxArray * fvec);
extern void mlxFdjac(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
