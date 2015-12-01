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

#include "libmatlb.h"
#include "fdjac.h"
#include "residual.h"

mxArray * Amp = NULL;

mxArray * C = NULL;

mxArray * M = NULL;

mxArray * Omega = NULL;

mxArray * P = NULL;

mxArray * Pamp = NULL;

mxArray * Pp = NULL;

mxArray * S = NULL;

mxArray * a = NULL;

mxArray * g = NULL;

mxArray * h0 = NULL;

mxArray * kappa = NULL;

mxArray * numunknowns = NULL;

mxArray * phi = NULL;

mxArray * w = NULL;

static mexGlobalTableEntry global_table[15]
  = { { "Amp", &Amp }, { "C", &C }, { "M", &M }, { "Omega", &Omega },
      { "P", &P }, { "Pamp", &Pamp }, { "Pp", &Pp }, { "S", &S },
      { "a", &a }, { "g", &g }, { "h0", &h0 }, { "kappa", &kappa },
      { "numunknowns", &numunknowns }, { "phi", &phi }, { "w", &w } };

static mexFunctionTableEntry function_table[2]
  = { { "fdjac", mlxFdjac, 2, 1, &_local_function_table_fdjac },
      { "residual", mlxResidual, 1, 1, NULL } };

static _mexInitTermTableEntry init_term_table[1]
  = { { InitializeModule_fdjac, TerminateModule_fdjac } };

static _mex_information _mex_info
  = { 1, 2, function_table, 15, global_table, 0, NULL, 1, init_term_table };

/*
 * The function "Mresidual" is the MATLAB callback version of the "residual"
 * function from file "C:\matlabR12\work\incompressible nonliear shallow water
 * equations\residual.m". It performs a callback to MATLAB to run the
 * "residual" function, and passes any resulting output arguments back to its
 * calling function.
 */
static mxArray * Mresidual(int nargout_, mxArray * coeffs) {
    mxArray * fvec = mclGetUninitializedArray();
    mclFevalCallMATLAB(
      mclNVarargout(nargout_, 0, &fvec, NULL), "residual", coeffs, NULL);
    return fvec;
}

/*
 * The function "mexLibrary" is a Compiler-generated mex wrapper, suitable for
 * building a MEX-function. It initializes any persistent variables as well as
 * a function table for use by the feval function. It then calls the function
 * "mlxFdjac". Finally, it clears the feval table and exits.
 */
mex_information mexLibrary(void) {
    return &_mex_info;
}

/*
 * The function "mlfResidual" contains the normal interface for the "residual"
 * M-function from file "C:\matlabR12\work\incompressible nonliear shallow
 * water equations\residual.m" (lines 0-0). This function processes any input
 * arguments and passes them to the implementation version of the function,
 * appearing above.
 */
mxArray * mlfResidual(mxArray * coeffs) {
    int nargout = 1;
    mxArray * fvec = mclGetUninitializedArray();
    mlfEnterNewContext(0, 1, coeffs);
    fvec = Mresidual(nargout, coeffs);
    mlfRestorePreviousContext(0, 1, coeffs);
    return mlfReturnValue(fvec);
}

/*
 * The function "mlxResidual" contains the feval interface for the "residual"
 * M-function from file "C:\matlabR12\work\incompressible nonliear shallow
 * water equations\residual.m" (lines 0-0). The feval function calls the
 * implementation version of residual through this function. This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
void mlxResidual(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: residual Line: 1 Column: "
            "1 The function \"residual\" was called with mor"
            "e than the declared number of outputs (1)."));
    }
    if (nrhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: residual Line: 1 Column:"
            " 1 The function \"residual\" was called with m"
            "ore than the declared number of inputs (1)."));
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = mclGetUninitializedArray();
    }
    for (i = 0; i < 1 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 1; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 1, mprhs[0]);
    mplhs[0] = Mresidual(nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
}
