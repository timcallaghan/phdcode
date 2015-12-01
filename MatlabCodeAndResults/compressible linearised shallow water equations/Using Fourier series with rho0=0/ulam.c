#include "mex.h"
#include <math.h>

/*
 * ulam.c
 *
 * Computational function that calculates the ulam
 * velocity component at a series of grid points in lat and long.
 *
 * This is a MEX-file for MATLAB.
 * 
 */


//*****************

double ulam1(double *P, int Plen, double phi)
{
	// Compute the order epsilon phi term
	double answer=0.0;
	int n,count;

	count = 0;

	for (n=1;n<=Plen;n++)
	{
		answer += *(P+count)*cos((2*n-1)*phi);
		count++;
	}

	return answer;

}

void ulam(double w, int kappa, double *eta, double *phi, double *P, int Plen, double *y, int mm, int nn, double epsilon)
{
	// Computational function that calculates the pressure
	// at a series of grid points in latitude and longitude.
	int ii,jj;

	for (ii=0; ii<mm; ii++) {
		for (jj=0; jj<nn; jj++) 
		{
			*(y+jj*mm+ii) = w*cos(*(phi+ii))+epsilon*cos(*(eta+jj)*kappa)*ulam1(P,Plen, *(phi+ii));
		}
	}
  
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *eta, *phi, *P, *y;
  int i,j,kappa,etarows,etacols,phirows,phicols,Prows,Pcols,mrows,ncols,Plen;
  double w,epsilon;
  
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
    mexErrMsgTxt("Input P must be a noncomplex double vector.");
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
    mexErrMsgTxt("Input epsilon must be a scalar.");
  }

  /*  get the scalar inputs w, a, kappa and epsilon*/
  kappa = mxGetScalar(prhs[3]);
  w = mxGetScalar(prhs[4]);
  epsilon = mxGetScalar(prhs[5]);

  /*  create a pointer to the input vectors eta, phi and D */
  eta = mxGetPr(prhs[0]);
  phi = mxGetPr(prhs[1]);
  P = mxGetPr(prhs[2]);

  /*  get the dimensions of the vectors eta, phi and D */
  etarows = mxGetM(prhs[0]);
  etacols = mxGetN(prhs[0]);
  phirows = mxGetM(prhs[1]);
  phicols = mxGetN(prhs[1]);
  Prows = mxGetM(prhs[2]);
  Pcols = mxGetN(prhs[2]);

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
  // Get the length of the vector D
  if (Prows > Pcols)
  {
	  Plen=Prows;
  }
  else
  {
	  Plen=Pcols;
  }
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
  
  /*  create a C pointer to a copy of the output matrix */
  y = mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  //pmnvec(n,m,phi,y,mrows,ncols);
  ulam(w,kappa,eta,phi,P,Plen,y,mrows,ncols,epsilon); 
}
