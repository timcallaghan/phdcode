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

double A(double w, double K, double Omega, double phi, int kappa)
{
	// Computes the A(phi) term of the free surface
	double temp1 = w*(2.0*Omega+w)*pow(cos(phi),2.0)/2.0;
	double temp2 = (K*K*pow(cos(phi),2.0*kappa)/4.0)*((kappa+1.0)*pow(cos(phi),2.0)+(2.0*kappa*kappa-kappa-2.0));
	double temp3 = -1.0*K*K*kappa*kappa*pow(cos(phi),2.0*kappa-2.0)/2.0;
	return temp1+temp2+temp3;
}

double B(double w, double K, double Omega, double phi, int kappa)
{
	// Computes the B(phi) term of the free surface
	double temp1 = (2.0*(Omega+w)*K/((kappa+1.0)*(kappa+2.0)))*pow(cos(phi),(double) kappa);
	double temp2 = (kappa*kappa+2.0*kappa+2.0)-(kappa+1.0)*(kappa+1.0)*pow(cos(phi),2.0);
	return temp1*temp2;
}

double C(double K, double phi, int kappa)
{
	// Computes the C(phi) term of the free surface
	double temp1 = K*K/4.0*pow(cos(phi),2.0*kappa);
	double temp2 = (kappa+1.0)*pow(cos(phi),2.0)-(kappa+2.0);
	return temp1*temp2;
}

void rhfsfield(double w, double a, double g, double Omega, int kappa, double *eta, double *phi, double *y, int mm, int nn, double K, double h0)
{
	// Computational function that calculates the pressure
	// at a series of grid points in latitude and longitude.
	int ii,jj;

	for (ii=0; ii<mm; ii++) {
		for (jj=0; jj<nn; jj++) 
		{
			*(y+jj*mm+ii) = h0 + a*a/g*A(w,K,Omega,*(phi+ii),kappa) + a*a/g*B(w,K,Omega,*(phi+ii),kappa)*cos( *(eta+jj)*kappa) + a*a/g*C(K,*(phi+ii),kappa)*cos( *(eta+jj)*kappa*2.0);
		}
	} 
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *eta, *phi, *y;
  int kappa,etarows,etacols,phirows,phicols,mrows,ncols;
  double w,a,g,Omega,K,h0;
  
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
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      mxGetN(prhs[3])*mxGetM(prhs[3])!=1 ) {
    mexErrMsgTxt("Input K must be a double scalar.");
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

  /* check to make sure the seventh input argument is a scalar */
  if( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) ||
      mxGetN(prhs[6])*mxGetM(prhs[6])!=1 ) {
    mexErrMsgTxt("Input g must be a scalar.");
  }

  /* check to make sure the eigth input argument is a scalar */
  if( !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) ||
      mxGetN(prhs[7])*mxGetM(prhs[7])!=1 ) {
    mexErrMsgTxt("Input Omega must be a scalar.");
  }

  /* check to make sure the ninth input argument is a scalar */
  if( !mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]) ||
      mxGetN(prhs[8])*mxGetM(prhs[8])!=1 ) {
    mexErrMsgTxt("Input h0 must be a scalar.");
  }

  /*  get the scalar inputs Fr, M, Ro, Sr, Ka, w, kappa and epsilon, h0*/
  K = mxGetScalar(prhs[2]);
  kappa = mxGetScalar(prhs[3]);
  w = mxGetScalar(prhs[4]);
  a = mxGetScalar(prhs[5]);
  g = mxGetScalar(prhs[6]);
  Omega = mxGetScalar(prhs[7]);
  h0 = mxGetScalar(prhs[8]);

  /*  create a pointer to the input vectors eta and phi */
  eta = mxGetPr(prhs[0]);
  phi = mxGetPr(prhs[1]);

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
  //pmnvec(n,m,phi,y,mrows,ncols);
  rhfsfield(w,a,g,Omega,kappa,eta,phi,y,mrows,ncols,K,h0); 
}
