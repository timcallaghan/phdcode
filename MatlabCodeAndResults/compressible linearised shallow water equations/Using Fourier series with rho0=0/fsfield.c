#include "mex.h"
#include <math.h>

/*
 * fsfield.c
 *
 * Computational function that calculates the free surface
 * field at a series of grid points in lat and long.
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


void fsfield(double w, double Fr, double Ro, int kappa, double *eta, double *phi, double *H, int Hlen, double *y, int mm, int nn, double epsilon, double h0)
{
	// Computational function that calculates the free surface height
	// at a series of grid points in latitude and longitude.
	int ii,jj;

	for (ii=0; ii<mm; ii++) {
		for (jj=0; jj<nn; jj++) 
		{
			*(y+jj*mm+ii) = hz(w,Fr,Ro, *(phi+ii),h0)+epsilon*cos( *(eta+jj)*kappa)*h1(H,Hlen, *(phi+ii));
		}
	}
  
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *eta, *phi, *H, *y;
  int i,j,kappa,etarows,etacols,phirows,phicols,Hrows,Hcols,mrows,ncols,Hlen;
  double w,Fr,Ro,epsilon,h0;
  
  /*  check for proper number of arguments */

  if(nrhs!=9) 
    mexErrMsgTxt("Nine inputs required.");
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
    mexErrMsgTxt("Input w must be a scalar.");
  }
  
  /* check to make sure the sixth input argument is a scalar */
  if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) ||
      mxGetN(prhs[5])*mxGetM(prhs[5])!=1 ) {
    mexErrMsgTxt("Input Fr must be a scalar.");
  }

  /* check to make sure the seventh input argument is a scalar */
  if( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) ||
      mxGetN(prhs[6])*mxGetM(prhs[6])!=1 ) {
    mexErrMsgTxt("Input Ro must be a scalar.");
  }

  /* check to make sure the eigth input argument is a scalar */
  if( !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) ||
      mxGetN(prhs[7])*mxGetM(prhs[7])!=1 ) {
    mexErrMsgTxt("Input epsilon must be a scalar.");
  }

  /* check to make sure the ninth input argument is a scalar */
  if( !mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]) ||
      mxGetN(prhs[8])*mxGetM(prhs[8])!=1 ) {
    mexErrMsgTxt("Input h0 must be a scalar.");
  }

  /*  get the scalar inputs kappa, w, Fr, Ro, epsilon and h0*/
  kappa = mxGetScalar(prhs[3]);
  w = mxGetScalar(prhs[4]);
  Fr = mxGetScalar(prhs[5]);
  Ro = mxGetScalar(prhs[6]);
  epsilon = mxGetScalar(prhs[7]);
  h0 = mxGetScalar(prhs[8]);

  /*  create a pointer to the input vectors eta, phi and H */
  eta = mxGetPr(prhs[0]);
  phi = mxGetPr(prhs[1]);
  H = mxGetPr(prhs[2]);

  /*  get the dimensions of the vectors eta, phi and H */
  etarows = mxGetM(prhs[0]);
  etacols = mxGetN(prhs[0]);
  phirows = mxGetM(prhs[1]);
  phicols = mxGetN(prhs[1]);
  Hrows = mxGetM(prhs[2]);
  Hcols = mxGetN(prhs[2]);

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
  if (Hrows > Hcols)
  {
	  Hlen=Hrows;
  }
  else
  {
	  Hlen=Hcols;
  }
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
  
  /*  create a C pointer to a copy of the output matrix */
  y = mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  fsfield(w,Fr,Ro,kappa,eta,phi,H,Hlen,y,mrows,ncols,epsilon,h0); 
}
