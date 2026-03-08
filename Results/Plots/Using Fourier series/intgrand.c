#include "mex.h"
#include <math.h>

/*
 * intgrand.c
 *
 * Computational function that calculates the integrand
 * of the linearized volume integration without simplification.
 *
 * This is a MEX-file for MATLAB.
 * 
 */

double h1(double *H, int Hlen, double phi)
{
	// Compute the order epsilon phi term
	double answer=0.0;
	int n,count;

	count = 0;

	for (n=1;n<=Hlen;n++)
	{
		answer += *(H+count)*(cos(2.0*n*phi)*pow(-1.0,(double) n)-cos(2.0*(n-1)*phi)*pow(-1.0,(double) n-1));
		count++;
	}

	return answer;

}

double hz(double w, double Fr, double Ro, double phi, double h0)
{
	// Computes the zonal flow free surface deviation
	
	double temp=w*Fr*Fr*(1.0/Ro+w)/2.0;
	return temp*pow(cos(phi),2.0)+h0;
}

void intgrand(double w, double Fr, double Ro, double h0, double a, double epsilon, int kappa, double *phi, double eta, double *H, int Hlen, double *y, int mm)
{
	// Computational function that calculates the integrand
	// of the linearized volume integration at a series of 
	// grid points in latitude phi.
	int ii;
	double tempval=0.0;

	for (ii=0; ii<mm; ii++) 
	{
		tempval=hz(w,Fr,Ro, *(phi+ii),h0)+epsilon*cos(kappa*eta)*h1(H,Hlen, *(phi+ii));
		*(y+ii) = (pow(tempval,3.0)+3.0*a*a*tempval+3.0*a*pow(tempval,2.0))*cos(*(phi+ii))/3.0;
	}
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *phi, *H, *y;
  int phirows,phicols,Hrows,Hcols,mrows,Hlen,kappa;
  double w,Fr,Ro,a,h0,epsilon,eta;
  
  /*  check for proper number of arguments */

  if(nrhs!=10) 
    mexErrMsgTxt("Ten inputs required.");
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
    mexErrMsgTxt("Input w must be a scalar.");
  }

  /* check to make sure the fifth input argument is a scalar */
  if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) ||
      mxGetN(prhs[4])*mxGetM(prhs[4])!=1 ) {
    mexErrMsgTxt("Input Fr must be a scalar.");
  }
  
  /* check to make sure the sixth input argument is a scalar */
  if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) ||
      mxGetN(prhs[5])*mxGetM(prhs[5])!=1 ) {
    mexErrMsgTxt("Input Ro must be a scalar.");
  }

  /* check to make sure the seventh input argument is a scalar */
  if( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) ||
      mxGetN(prhs[6])*mxGetM(prhs[6])!=1 ) {
    mexErrMsgTxt("Input a must be a scalar.");
  }
  
  /* check to make sure the eighth input argument is a scalar */
  if( !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) ||
      mxGetN(prhs[7])*mxGetM(prhs[7])!=1 ) {
    mexErrMsgTxt("Input h0 must be a scalar.");
  }

  /* check to make sure the nineth input argument is a scalar */
  if( !mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]) ||
      mxGetN(prhs[8])*mxGetM(prhs[8])!=1 ) {
    mexErrMsgTxt("Input kappa must be a scalar (integer).");
  }

  /* check to make sure the tenth input argument is a scalar */
  if( !mxIsDouble(prhs[9]) || mxIsComplex(prhs[9]) ||
      mxGetN(prhs[9])*mxGetM(prhs[9])!=1 ) {
    mexErrMsgTxt("Input epsilon must be a scalar.");
  }

  /*  get the scalar inputs kappa, w, Fr, Ro, epsilon and h0*/
  eta = mxGetScalar(prhs[1]);
  w = mxGetScalar(prhs[3]);
  Fr = mxGetScalar(prhs[4]);
  Ro = mxGetScalar(prhs[5]);
  a = mxGetScalar(prhs[6]);
  h0 = mxGetScalar(prhs[7]);
  kappa = mxGetScalar(prhs[8]);
  epsilon = mxGetScalar(prhs[9]);
  

  /*  create a pointer to the input vectors phi and H */
  phi = mxGetPr(prhs[0]);
  H = mxGetPr(prhs[2]);

  /*  get the dimensions of the vectors eta, phi and H */
  phirows = mxGetM(prhs[0]);
  phicols = mxGetN(prhs[0]);
  Hrows = mxGetM(prhs[2]);
  Hcols = mxGetN(prhs[2]);

  // Get the length of the vector phi
  if (phirows > phicols)
  {
	  mrows=phirows;
  }
  else
  {
	  mrows=phicols;
  }
  // Get the length of the vector H
  if (Hrows > Hcols)
  {
	  Hlen=Hrows;
  }
  else
  {
	  Hlen=Hcols;
  }
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(mrows,1, mxREAL);
  
  /*  create a C pointer to a copy of the output matrix */
  y = mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  intgrand(w,Fr,Ro,h0,a,epsilon,kappa,phi,eta,H,Hlen,y,mrows); 
}