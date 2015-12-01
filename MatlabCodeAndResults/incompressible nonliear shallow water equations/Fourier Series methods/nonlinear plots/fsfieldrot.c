#include "mex.h"
#include <math.h>

/*
 * fsfieldrot.c
 *
 * Computational function that calculates the nonlinear free surface
 * field at a series of grid points in lat and long and rotates it by a
 * certain angle
 *
 * This is a MEX-file for MATLAB.
 * 
 */

double h1(double *H, double phi, double eta, int M, int N, int kappa, double rotangle)
{
	// Compute the nonlinear series at phi and eta
	double answer=0.0;
	int n,m,count;
	// First compute the eta-independent terms
	for (n=0;n<=N;n++)
	{
		answer += *(H+n)*cos(2.0*n*phi);
	}
	// Now compute the rest of the series
	count=N+1;
	for (m=1;m<=M-1;m++)
	{
		for (n=1;n<=N;n++)
		{
			answer += *(H+count)*cos(kappa*m*(eta-rotangle))*pow(-1.0,(double) n)*(cos(2.0*n*phi)+cos(2.0*(n-1)*phi));
			count++;
		}
	}

	return answer;
}

void fsfield(int M, int N, int kappa, double *eta, double *phi, double *H, double *y, int mm, int nn, double rotangle)
{
	// Computational function that calculates the free surface height
	// at a series of grid points in latitude and longitude.
	int ii,jj;

	for (ii=0; ii<mm; ii++) {
		for (jj=0; jj<nn; jj++) 
		{
			*(y+jj*mm+ii) = h1(H,*(phi+ii),*(eta+jj),M,N,kappa,rotangle);
		}
	}
  
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *eta, *phi, *H, *y;
  int M,N,kappa,etarows,etacols,phirows,phicols,mrows,ncols;
  double rotangle;
  
  /*  check for proper number of arguments */
  if(nrhs!=7) 
    mexErrMsgTxt("Seven inputs required.");
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
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ) {
    mexErrMsgTxt("Input H must be a noncomplex double vector.");
  }

  /* check to make sure the fourth input argument is a scalar */
  if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
      mxGetN(prhs[3])*mxGetM(prhs[3])!=1 ) {
    mexErrMsgTxt("Input kappa must be a scalar(integer).");
  }

  /* check to make sure the fifth input argument is a scalar */
  if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) ||
      mxGetN(prhs[4])*mxGetM(prhs[4])!=1 ) {
    mexErrMsgTxt("Input M must be a scalar(integer).");
  }
  
  /* check to make sure the sixth input argument is a scalar */
  if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) ||
      mxGetN(prhs[5])*mxGetM(prhs[5])!=1 ) {
    mexErrMsgTxt("Input N must be a scalar(integer).");
  }

    /* check to make sure the seventh input argument is a scalar */
  if( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) ||
      mxGetN(prhs[6])*mxGetM(prhs[6])!=1 ) {
    mexErrMsgTxt("Input rotangle must be a scalar.");
  }

  /*  get the scalar inputs kappa, M and N*/
  kappa = mxGetScalar(prhs[3]);
  M = mxGetScalar(prhs[4]);
  N = mxGetScalar(prhs[5]);
  rotangle = mxGetScalar(prhs[6]);

  /*  create a pointer to the input vectors eta, phi and H */
  eta = mxGetPr(prhs[0]);
  phi = mxGetPr(prhs[1]);
  H = mxGetPr(prhs[2]);

  /*  get the dimensions of the vectors eta and phi */
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
  fsfield(M,N,kappa,eta,phi,H,y,mrows,ncols,rotangle); 
}
