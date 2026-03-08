#include "mex.h"
#include <math.h>

/*
 * I3.c
 *
 * Computational function that calculates the integrand
 * of the derivative of the nonlinear volume integration.
 * with repsect to a coefficient.
 *
 * This is a MEX-file for MATLAB.
 * 
 */

double h1(double *H, double phi,double eta, int M, int N, int kappa)
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
			answer += *(H+count)*cos(kappa*m*eta)*pow(-1.0,(double) n)*(cos(2.0*n*phi)+cos(2.0*(n-1)*phi));
			count++;
		}
	}

	return answer;
}

void I3(int M, int N, double a, int kappa, double *phi, double eta, double *H, double *y, int mm, int nval, int mval)
{
	// Computational function that calculates the integrand
	// of the nonlinear volume integration at a series of 
	// grid points in latitude phi.
	int ii;
	double tempval=0.0;
	double dhdHval = 0.0;
	double phival = 0.0;

	for (ii=0; ii<mm; ii++) 
	{
		phival = *(phi+ii);
		tempval=h1(H,phival,eta,M,N,kappa);
		if (mval==0)
		{
			dhdHval=cos(phival*2.0*nval);
		}
		else
		{
			dhdHval=cos(kappa*mval*eta)*pow(-1.0,(double) nval)*(cos(2.0*nval*phival)+cos(2.0*(nval-1)*phival));
		}
		*(y+ii) = pow(tempval,2.0)*cos(phival)*dhdHval;
	}
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *phi, *H, *y;
  int M,N,phirows,phicols,mrows,kappa,nval,mval;
  double a,eta;
  
  /*  check for proper number of arguments */

  if(nrhs!=9) 
    mexErrMsgTxt("Seven inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");

  /* check to make sure the first input argument is a noncomplex
  double (can be scalar or vector*/
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) {
    mexErrMsgTxt("Input phi must be a noncomplex double vector.");
  }

  /* check to make sure the fourth input argument is a scalar */
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
      mxGetN(prhs[1])*mxGetM(prhs[1])!=1 ) {
    mexErrMsgTxt("Input eta must be a scalar.");
  }

  /* check to make sure the second input argument is a noncomplex
  double (can be scalar or vector*/
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ) {
    mexErrMsgTxt("Input H must be a noncomplex double vector.");
  }

  /* check to make sure the fourth input argument is a scalar */
  if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
      mxGetN(prhs[3])*mxGetM(prhs[3])!=1 ) {
    mexErrMsgTxt("Input M must be a scalar(integer).");
  }

  /* check to make sure the fifth input argument is a scalar */
  if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) ||
      mxGetN(prhs[4])*mxGetM(prhs[4])!=1 ) {
    mexErrMsgTxt("Input N must be a scalar(integer).");
  }
  
  /* check to make sure the sixth input argument is a scalar */
  if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) ||
      mxGetN(prhs[5])*mxGetM(prhs[5])!=1 ) {
    mexErrMsgTxt("Input kappa must be a scalar(integer).");
  }

  /* check to make sure the seventh input argument is a scalar */
  if( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) ||
      mxGetN(prhs[6])*mxGetM(prhs[6])!=1 ) {
    mexErrMsgTxt("Input a must be a scalar.");
  }

  /* check to make sure the eigth input argument is a scalar */
  if( !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) ||
      mxGetN(prhs[7])*mxGetM(prhs[7])!=1 ) {
    mexErrMsgTxt("Input nval must be a scalar(integer).");
  }  

  /* check to make sure the ninth input argument is a scalar */
  if( !mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]) ||
      mxGetN(prhs[8])*mxGetM(prhs[8])!=1 ) {
    mexErrMsgTxt("Input mval must be a scalar(integer).");
  } 
  
  /*  get the scalar inputs kappa, w, Fr, Ro, epsilon and h0*/
  eta = mxGetScalar(prhs[1]);
  M = mxGetScalar(prhs[3]);
  N = mxGetScalar(prhs[4]);
  kappa = mxGetScalar(prhs[5]);
  a = mxGetScalar(prhs[6]);
  nval = mxGetScalar(prhs[7]);
  mval = mxGetScalar(prhs[8]);

  /*  create a pointer to the input vectors phi and H */
  phi = mxGetPr(prhs[0]);
  H = mxGetPr(prhs[2]);

  /*  get the dimensions of the vectors phi and H */
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
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(mrows,1, mxREAL);
  
  /*  create a C pointer to a copy of the output matrix */
  y = mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  I3(M,N,a,kappa,phi,eta,H,y,mrows,nval,mval); 
}