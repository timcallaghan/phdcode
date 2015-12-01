/*
 * MATLAB Compiler: 2.1
 * Date: Tue Mar 09 12:17:15 2004
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C"
 * "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "fdjac" 
 */
#include "fdjac.h"
#include "libmatlbm.h"
#include "residual.h"

extern mxArray * numunknowns;

static mxChar _array1_[128] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 'f', 'd', 'j', 'a', 'c',
                                ' ', 'L', 'i', 'n', 'e', ':', ' ', '1', ' ',
                                'C', 'o', 'l', 'u', 'm', 'n', ':', ' ', '1',
                                ' ', 'T', 'h', 'e', ' ', 'f', 'u', 'n', 'c',
                                't', 'i', 'o', 'n', ' ', '"', 'f', 'd', 'j',
                                'a', 'c', '"', ' ', 'w', 'a', 's', ' ', 'c',
                                'a', 'l', 'l', 'e', 'd', ' ', 'w', 'i', 't',
                                'h', ' ', 'm', 'o', 'r', 'e', ' ', 't', 'h',
                                'a', 'n', ' ', 't', 'h', 'e', ' ', 'd', 'e',
                                'c', 'l', 'a', 'r', 'e', 'd', ' ', 'n', 'u',
                                'm', 'b', 'e', 'r', ' ', 'o', 'f', ' ', 'o',
                                'u', 't', 'p', 'u', 't', 's', ' ', '(', '1',
                                ')', '.' };
static mxArray * _mxarray0_;

static mxChar _array3_[127] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 'f', 'd', 'j', 'a', 'c',
                                ' ', 'L', 'i', 'n', 'e', ':', ' ', '1', ' ',
                                'C', 'o', 'l', 'u', 'm', 'n', ':', ' ', '1',
                                ' ', 'T', 'h', 'e', ' ', 'f', 'u', 'n', 'c',
                                't', 'i', 'o', 'n', ' ', '"', 'f', 'd', 'j',
                                'a', 'c', '"', ' ', 'w', 'a', 's', ' ', 'c',
                                'a', 'l', 'l', 'e', 'd', ' ', 'w', 'i', 't',
                                'h', ' ', 'm', 'o', 'r', 'e', ' ', 't', 'h',
                                'a', 'n', ' ', 't', 'h', 'e', ' ', 'd', 'e',
                                'c', 'l', 'a', 'r', 'e', 'd', ' ', 'n', 'u',
                                'm', 'b', 'e', 'r', ' ', 'o', 'f', ' ', 'i',
                                'n', 'p', 'u', 't', 's', ' ', '(', '2', ')',
                                '.' };
static mxArray * _mxarray2_;
static mxArray * _mxarray4_;
static mxArray * _mxarray5_;
static mxArray * _mxarray6_;
static mxArray * _mxarray7_;

void InitializeModule_fdjac(void) {
    _mxarray0_ = mclInitializeString(128, _array1_);
    _mxarray2_ = mclInitializeString(127, _array3_);
    _mxarray4_ = mclInitializeDouble(1e-08);
    _mxarray5_ = mclInitializeDouble(0.0);
    _mxarray6_ = mclInitializeDouble(1.0);
    _mxarray7_ = mclInitializeDoubleVector(0, 0, (double *)NULL);
}

void TerminateModule_fdjac(void) {
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mfdjac(int nargout_, mxArray * x, mxArray * fvec);

_mexLocalFunctionTable _local_function_table_fdjac
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfFdjac" contains the normal interface for the "fdjac"
 * M-function from file "C:\matlabR12\work\incompressible nonliear shallow
 * water equations\fdjac.m" (lines 1-34). This function processes any input
 * arguments and passes them to the implementation version of the function,
 * appearing above.
 */
mxArray * mlfFdjac(mxArray * x, mxArray * fvec) {
    int nargout = 1;
    mxArray * df = mclGetUninitializedArray();
    mlfEnterNewContext(0, 2, x, fvec);
    df = Mfdjac(nargout, x, fvec);
    mlfRestorePreviousContext(0, 2, x, fvec);
    return mlfReturnValue(df);
}

/*
 * The function "mlxFdjac" contains the feval interface for the "fdjac"
 * M-function from file "C:\matlabR12\work\incompressible nonliear shallow
 * water equations\fdjac.m" (lines 1-34). The feval function calls the
 * implementation version of fdjac through this function. This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
void mlxFdjac(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(_mxarray0_);
    }
    if (nrhs > 2) {
        mlfError(_mxarray2_);
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = mclGetUninitializedArray();
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    mplhs[0] = Mfdjac(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mfdjac" is the implementation version of the "fdjac"
 * M-function from file "C:\matlabR12\work\incompressible nonliear shallow
 * water equations\fdjac.m" (lines 1-34). It contains the actual compiled code
 * for that M-function. It is a static function and must only be called from
 * one of the interface functions, appearing below.
 */
/*
 * function df = fdjac(x,fvec)
 */
static mxArray * Mfdjac(int nargout_, mxArray * x, mxArray * fvec) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_fdjac);
    mxArray * df = mclGetUninitializedArray();
    mxArray * j = mclGetUninitializedArray();
    mxArray * fd = mclGetUninitializedArray();
    mxArray * f = mclGetUninitializedArray();
    mxArray * n = mclGetUninitializedArray();
    mxArray * temp = mclGetUninitializedArray();
    mxArray * h = mclGetUninitializedArray();
    mxArray * EPS = mclGetUninitializedArray();
    mxArray * ans = mclGetUninitializedArray();
    mclCopyArray(&x);
    mclCopyArray(&fvec);
    /*
     * % Computes forward-difference approximation to Jacobian. On
     * % input, x[1..n] is the point at which the Jacobian is to be
     * % evaluated, fvec[1..n] is the vector of function values at
     * % the point, and residual(x) is a user-supplied routine that
     * % returns the vector of functions at x. On output, df[1..n][1..n]
     * % is the Jaconian array.
     * 
     * % Global variables
     * global M kappa g a Omega w h0 phi P Pp C S Pamp numunknowns Amp
     * 
     * % Approximate square root of the machine precision.
     * EPS=1.0e-8; 
     */
    mlfAssign(&EPS, _mxarray4_);
    /*
     * h=0;
     */
    mlfAssign(&h, _mxarray5_);
    /*
     * temp=0;
     */
    mlfAssign(&temp, _mxarray5_);
    /*
     * n=numunknowns;
     */
    mlfAssign(&n, mclVg(&numunknowns, "numunknowns"));
    /*
     * f=zeros(n,1);
     */
    mlfAssign(&f, mlfZeros(mclVv(n, "n"), _mxarray6_, NULL));
    /*
     * fd=zeros(n,n);
     */
    mlfAssign(&fd, mlfZeros(mclVv(n, "n"), mclVv(n, "n"), NULL));
    /*
     * for j=1:n
     */
    {
        int v_ = mclForIntStart(1);
        int e_ = mclForIntEnd(mclVv(n, "n"));
        if (v_ > e_) {
            mlfAssign(&j, _mxarray7_);
        } else {
            /*
             * temp=x(j,1);
             * h=EPS*abs(temp);
             * if (h == 0.0) 
             * h=EPS;
             * end
             * % Trick to reduce finite precision error.
             * x(j,1)=temp+h;
             * h=x(j,1)-temp;
             * f=residual(x);
             * x(j,1)=temp;
             * % Calculate the FD jacobian column-wise
             * df(:,j)=(f-fvec)/h;
             * end
             */
            for (; ; ) {
                mlfAssign(&temp, mclIntArrayRef2(mclVsa(x, "x"), v_, 1));
                mlfAssign(
                  &h,
                  mclMtimes(
                    mclVv(EPS, "EPS"), mclVe(mlfAbs(mclVv(temp, "temp")))));
                if (mclEqBool(mclVv(h, "h"), _mxarray5_)) {
                    mlfAssign(&h, mclVsv(EPS, "EPS"));
                }
                mclIntArrayAssign2(
                  &x, mclPlus(mclVv(temp, "temp"), mclVv(h, "h")), v_, 1);
                mlfAssign(
                  &h,
                  mclMinus(
                    mclVe(mclIntArrayRef2(mclVsa(x, "x"), v_, 1)),
                    mclVv(temp, "temp")));
                mlfAssign(&f, mlfResidual(mclVa(x, "x")));
                mclIntArrayAssign2(&x, mclVsv(temp, "temp"), v_, 1);
                mclArrayAssign2(
                  &df,
                  mclMrdivide(
                    mclMinus(mclVv(f, "f"), mclVa(fvec, "fvec")),
                    mclVv(h, "h")),
                  mlfCreateColonIndex(),
                  mlfScalar(v_));
                if (v_ == e_) {
                    break;
                }
                ++v_;
            }
            mlfAssign(&j, mlfScalar(v_));
        }
    }
    mclValidateOutput(df, 1, nargout_, "df", "fdjac");
    mxDestroyArray(ans);
    mxDestroyArray(EPS);
    mxDestroyArray(h);
    mxDestroyArray(temp);
    mxDestroyArray(n);
    mxDestroyArray(f);
    mxDestroyArray(fd);
    mxDestroyArray(j);
    mxDestroyArray(fvec);
    mxDestroyArray(x);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return df;
    /*
     * 
     */
}
