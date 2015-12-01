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

#ifndef __residual_h
#define __residual_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_residual(void);
extern void TerminateModule_residual(void);
extern _mexLocalFunctionTable _local_function_table_residual;

extern mxArray * mlfResidual(mxArray * coeffs);
extern void mlxResidual(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
