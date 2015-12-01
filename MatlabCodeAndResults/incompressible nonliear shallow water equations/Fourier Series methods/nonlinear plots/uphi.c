#include "mex.h"
#include <math.h>

/*
 * uphi.c
 *
 * Computational function that calculates the uphi
 * velocity component at a series of grid points in lat and long.
 *
 * This is a MEX-file for MATLAB.
 * 
 */
double uphi1(double *Q, double phi, double eta, int M, int N, int kappa)
{
	// Compute the nonlinear series at phi and eta
	double answer=0.0;
	int n,m,count;
	// Now compute the rest of the series
	count=0;
	for (m=1;m<=M;m++)
	{
		for (n=1;n<=N;n++)
		{
			answer += *(Q+count)*sin(kappa*m*eta)*sin(2.0*n*phi);
			count++;
		}
	}

	return answer;
}

void uphi(int M, int N, int kappa, double *eta, double *phi, double *Q, double *y, int mm, int nn)
{
	// Computational function that calculates the pressure
	// at a series of grid points in latitude and longitude.
	int ii,jj;

	for (ii=0; ii<mm; ii++) {
		for (jj=0; jj<nn; jj++) 
		{
			*(y+jj*mm+ii) = uphi1(Q,*(phi+ii),*(eta+jj),M,N,kappa);
		}
	}
  
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *eta, *phi, *Q, *y;
  int M,N,kappa,etarows,etacols,phirows,phicols,mrows,ncols;
  
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
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ) {
    mexErrMsgTxt("Input Q must be a noncomplex double vector.");
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

  /* check to make sure the fifth input argument is a scalar */
  if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) ||
      mxGetN(prhs[5])*mxGetM(prhs[5])!=1 ) {
    mexErrMsgTxt("Input N must be a scalar(integer).");
  }

  /*  get the scalar inputs kappa and epsilon*/
  kappa = mxGetScalar(prhs[3]);
  M = mxGetScalar(prhs[4]);
  N = mxGetScalar(prhs[5]);

  /*  create a pointer to the input vectors eta, phi and B */
  eta = mxGetPr(prhs[0]);
  phi = mxGetPr(prhs[1]);
  Q = mxGetPr(prhs[2]);

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
  uphi(M,N,kappa,eta,phi,Q,y,mrows,ncols); 
}
