#include "mex.h"
#include <math.h>

/*
 * rhulam.c
 *
 * Computational function that calculates the ulam
 * velocity component at a series of grid points in lat and long.
 *
 * This is a MEX-file for MATLAB.
 * 
 */

//*****************


void rhulam(double K, double w, double a, int kappa, double *eta, double *phi, double *y, int mm, int nn)
{
	// Computational function that calculates the pressure
	// at a series of grid points in latitude and longitude.
	int ii,jj;

	for (ii=0; ii<mm; ii++) {
		for (jj=0; jj<nn; jj++) 
		{
			*(y+jj*mm+ii) = w*a*cos(*(phi+ii))+a*K*pow(cos(*(phi+ii)),kappa-1.0)*(kappa*pow(sin(*(phi+ii)),2.0)-pow(cos(*(phi+ii)),2.0))*cos(*(eta+jj)*kappa);
		}
	}
  
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *eta, *phi, *y;
  int i,j,kappa,etarows,etacols,phirows,phicols,mrows,ncols;
  double w,a,K;
  
  /*  check for proper number of arguments */

  if(nrhs!=6) 
    mexErrMsgTxt("Six inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");

  /* check to make sure the first input argument is a noncomplex
  double (can be scalar or vector*/
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) {
    mexErrMsgTxt("Input eta must be a noncomplex double vector.");
  }

  /* check to make sure the second input argument is a noncomplex
  double (can be scalar or vector*/
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) {
    mexErrMsgTxt("Input phi must be a noncomplex double vector.");
  }

  /* check to make sure the third input argument is a noncomplex
  double (can be scalar or vector*/
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      mxGetN(prhs[3])*mxGetM(prhs[3])!=1 ) {
    mexErrMsgTxt("Input K must be a noncomplex double vector.");
  }

  /* check to make sure the fourth input argument is a scalar */
  if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
      mxGetN(prhs[3])*mxGetM(prhs[3])!=1 ) {
    mexErrMsgTxt("Input kappa must be a scalar(integer).");
  }

  /* check to make sure the fifth input argument is a scalar */
  if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) ||
      mxGetN(prhs[4])*mxGetM(prhs[4])!=1 ) {
    mexErrMsgTxt("Input w must be a scalar.");
  }
  
  /* check to make sure the sixth input argument is a scalar */
  if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) ||
      mxGetN(prhs[5])*mxGetM(prhs[5])!=1 ) {
    mexErrMsgTxt("Input a must be a scalar.");
  }

  /*  get the scalar inputs w, a, kappa and K*/
  K = mxGetScalar(prhs[2]);
  kappa = mxGetScalar(prhs[3]);
  w = mxGetScalar(prhs[4]);
  a = mxGetScalar(prhs[5]);

  /*  create a pointer to the input vectors eta and phi*/
  eta = mxGetPr(prhs[0]);
  phi = mxGetPr(prhs[1]);

  /*  get the dimensions of the vectors eta, phi and D */
  etarows = mxGetM(prhs[0]);
  etacols = mxGetN(prhs[0]);
  phirows = mxGetM(prhs[1]);
  phicols = mxGetN(prhs[1]);

  // Get the length of the vector eta
  if (etarows > etacols)
  {
	  ncols=etarows;
  }
  else
  {
	  ncols=etacols;
  }
  // Get the length of the vector phi
  if (phirows > phicols)
  {
	  mrows=phirows;
  }
  else
  {
	  mrows=phicols;
  }
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
  
  /*  create a C pointer to a copy of the output matrix */
  y = mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  //pmnvec(n,m,phi,y,mrows,ncols);
  rhulam(K,w,a,kappa,eta,phi,y,mrows,ncols); 
}
