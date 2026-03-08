#include "mex.h"
#include <math.h>

/*
 * integrnd.c
 *
 * Computational function that calculates the integrand
 * of the linearized volume integration.
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

void integrnd(double w, double Fr, double Ro, double h0, double a, double *phi, double *H, int Hlen, double *y, int mm)
{
	// Computational function that calculates the integrand
	// of the linearized volume integration at a series of 
	// grid points in latitude phi.
	int ii;
	double tempval=0.0;

	for (ii=0; ii<mm; ii++) 
	{
		tempval=h1(H,Hlen, *(phi+ii));
		*(y+ii) = (hz(w,Fr,Ro, *(phi+ii),h0)+a)*cos(*(phi+ii))*pow(tempval,2.0);
	}
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *phi, *H, *y;
  int phirows,phicols,Hrows,Hcols,mrows,Hlen;
  double w,Fr,Ro,a,h0;
  
  /*  check for proper number of arguments */

  if(nrhs!=7) 
    mexErrMsgTxt("Seven inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");

  /* check to make sure the first input argument is a noncomplex
  double (can be scalar or vector*/
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) {
    mexErrMsgTxt("Input phi must be a noncomplex double vector.");
  }

  /* check to make sure the second input argument is a noncomplex
  double (can be scalar or vector*/
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) {
    mexErrMsgTxt("Input H must be a noncomplex double vector.");
  }

  /* check to make sure the third input argument is a scalar */
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      mxGetN(prhs[2])*mxGetM(prhs[2])!=1 ) {
    mexErrMsgTxt("Input h0 must be a scalar double.");
  }

  /* check to make sure the fourth input argument is a scalar */
  if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
      mxGetN(prhs[3])*mxGetM(prhs[3])!=1 ) {
    mexErrMsgTxt("Input w must be a scalar(integer).");
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

  /*  get the scalar inputs kappa, w, Fr, Ro, epsilon and h0*/
  h0 = mxGetScalar(prhs[2]);
  w = mxGetScalar(prhs[3]);
  Fr = mxGetScalar(prhs[4]);
  Ro = mxGetScalar(prhs[5]);
  a = mxGetScalar(prhs[6]);
  

  /*  create a pointer to the input vectors phi and H */
  phi = mxGetPr(prhs[0]);
  H = mxGetPr(prhs[1]);

  /*  get the dimensions of the vectors eta, phi and H */
  phirows = mxGetM(prhs[0]);
  phicols = mxGetN(prhs[0]);
  Hrows = mxGetM(prhs[1]);
  Hcols = mxGetN(prhs[1]);

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
  integrnd(w,Fr,Ro,h0,a,phi,H,Hlen,y,mrows); 
}

