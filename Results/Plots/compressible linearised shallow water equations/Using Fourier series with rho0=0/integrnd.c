#include "mex.h"
#include <math.h>

/*
 * integrnd.c
 *
 * Computational function that calculates the integrand
 * of the linearized compressible mass integration.
 *
 * This is a MEX-file for MATLAB.
 * 
 */

void integrnd(double h0, double beta, double Aw, double a, double gamma, double *phi, double *y, int mm)
{
	// Computational function that calculates the integrand
	// of the linearized compressible mass integration at a series of 
	// grid points in latitude phi.
	int ii;
	double tempval=0.0;
	double f1val = 0.0;
	double hzval = 0.0;

	for (ii=0; ii<mm; ii++) 
	{
		hzval = h0 + Aw*pow(cos(*(phi+ii)),2.0);
		f1val = beta*(a+hzval);
		tempval = pow(hzval,gamma/(gamma-1.0));
		*(y+ii) = cos(*(phi+ii))*tempval*(2.0*pow(f1val*(gamma-1.0)/a,2.0)+2.0*beta*(gamma-1.0)*gamma*(f1val/a)+gamma*(2.0*gamma-1.0)*pow(beta,2.0));
	}
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *phi, *y;
  int phirows,phicols,mrows;
  double h0,beta,Aw,a,gamma;
  
  /*  check for proper number of arguments */

  if(nrhs!=6) 
    mexErrMsgTxt("Six inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");

  /* check to make sure the first input argument is a noncomplex
  double (can be scalar or vector*/
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) {
    mexErrMsgTxt("Input phi must be a noncomplex double vector.");
  }

  /* check to make sure the second input argument is a scalar */
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
      mxGetN(prhs[1])*mxGetM(prhs[1])!=1 ) {
    mexErrMsgTxt("Input h0 must be a scalar double.");
  }

  /* check to make sure the third input argument is a scalar */
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      mxGetN(prhs[2])*mxGetM(prhs[2])!=1 ) {
    mexErrMsgTxt("Input beta must be a scalar double.");
  }

  /* check to make sure the fourth input argument is a scalar */
  if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
      mxGetN(prhs[3])*mxGetM(prhs[3])!=1 ) {
    mexErrMsgTxt("Input Aw must be a scalar double.");
  }
  
  /* check to make sure the fifth input argument is a scalar */
  if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) ||
      mxGetN(prhs[4])*mxGetM(prhs[4])!=1 ) {
    mexErrMsgTxt("Input a must be a scalar.");
  }

  /* check to make sure the sixth input argument is a scalar */
  if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) ||
      mxGetN(prhs[5])*mxGetM(prhs[5])!=1 ) {
    mexErrMsgTxt("Input gamma must be a scalar.");
  }

  /*  get the scalar inputs h0, beta, Aw, a and gamma*/
  h0 = mxGetScalar(prhs[1]);
  beta = mxGetScalar(prhs[2]);
  Aw = mxGetScalar(prhs[3]);
  a = mxGetScalar(prhs[4]);
  gamma = mxGetScalar(prhs[5]);
  

  /*  create a pointer to the input vector phi */
  phi = mxGetPr(prhs[0]);

  /*  get the dimensions of the vector phi */
  phirows = mxGetM(prhs[0]);
  phicols = mxGetN(prhs[0]);

  // Get the length of the vector phi
  if (phirows > phicols)
  {
	  mrows=phirows;
  }
  else
  {
	  mrows=phicols;
  }
  
  /*  set the output pointer to the output vector (only 1 column) */
  plhs[0] = mxCreateDoubleMatrix(mrows,1, mxREAL);
  
  /*  create a C pointer to a copy of the output matrix */
  y = mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  integrnd(h0,beta,Aw,a,gamma,phi,y,mrows); 
}

