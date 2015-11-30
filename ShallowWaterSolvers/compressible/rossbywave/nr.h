#ifndef _NR_H_
#define _NR_H_
#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include "nrutil.h"
#include "nrtypes.h"
using namespace std;

namespace NR {
void cacheCeta(Mat_O_DP &Ceta);
void cacheSeta(Mat_O_DP &Seta);
void cacheCphi1(Mat_O_DP &Cphi1, Vec_I_DP &phi);
void cacheSphi1(Mat_O_DP &Sphi1, Vec_I_DP &phi);
void cacheCphi2(Mat_O_DP &Cphi2, Vec_I_DP &phi);
void cacheSphi2(Mat_O_DP &Sphi2, Vec_I_DP &phi);
void mnewt(Vec_IO_DP &coeffs,Vec_I_DP &phi,Mat_I_DP &Ceta,Mat_I_DP &Seta,Mat_I_DP &Cphi1,Mat_I_DP &Sphi1,Mat_I_DP &Cphi2,Mat_I_DP &Sphi2);
void residual(Vec_IO_DP &coeffs,Vec_I_DP &phi,Mat_I_DP &Ceta,Mat_I_DP &Seta,Mat_I_DP &Cphi1,Mat_I_DP &Sphi1,Mat_I_DP &Cphi2, Mat_I_DP &Sphi2,Vec_IO_DP &fvec);
void jacobian(Vec_IO_DP &coeffs,Vec_I_DP &phi,Mat_I_DP &Ceta,Mat_I_DP &Seta,Mat_I_DP &Cphi1,Mat_I_DP &Sphi1,Mat_I_DP &Cphi2, Mat_I_DP &Sphi2,Mat_IO_DP &jac);
//void fdjac(Vec_IO_DP &coeffs,Vec_I_DP &phi,Mat3D_I_DP &P,Mat3D_I_DP &Pp,Mat_I_DP &C,Mat_I_DP &S,Vec_O_DP &fvec,Mat_O_DP &df);
//void lubksb(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b);
//void ludcmp(Mat_IO_DP &a, Vec_O_INT &indx, DP &d);
DP quad2d(DP func(const DP, const DP, Vec_I_DP &coeffs), const DP x1, const DP x2, const DP yy1, const DP yy2, Vec_I_DP &coeffs, const DP tol);
DP quad2dmod(DP func(const DP, const DP, Vec_I_DP &coeffs, const int, const int), const DP x1, const DP x2, const DP yy1, const DP yy2, Vec_I_DP &coeffs, const DP tol, const int mval, const int nval);
DP adaptlob(DP f(const DP),DP const a,DP const b,DP tol);
void qrdcmp(Mat_IO_DP &a, Vec_O_DP &c, Vec_O_DP &d, bool &sing);
void qrsolv(Mat_I_DP &a, Vec_I_DP &c, Vec_I_DP &d, Vec_IO_DP &b);
void rsolv(Mat_I_DP &a, Vec_I_DP &d, Vec_IO_DP &b);
void loadLinCoeffs(Vec_IO_DP & coeffs);
void loadNLinCoeffs(Vec_IO_DP & coeffs);
}
#endif /* _NR_H_ */

