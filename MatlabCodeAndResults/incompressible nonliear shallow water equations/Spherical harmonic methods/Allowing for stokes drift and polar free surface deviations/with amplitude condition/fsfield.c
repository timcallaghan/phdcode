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



// Used in evaluating factorial n for integers n>32
double gammln(const double xx)
{
	// Returns the value of ln[gamma(xx)] for xx > 0;
	int j;
	double x,y,tmp,ser;
	static const double cof[6]={76.18009172947146, -86.50532032941677,
		24.01409824083091,-1.231739572450155,0.128650973866179e-2,
		-0.5395239384953e-5};

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<6;j++)
	{
		ser += cof[j]/++y;
	}
	return -tmp+log(2.5066282746310005*ser/x);
}

// Used to calculate factorial n for n<32
double factrl(const int n)
{
	// Returns the value of n! as a floating point number.
	static int ntop=4;
	// Fill in table only as required...
	static double a[33]={1.0,1.0,2.0,6.0,24.0};
	int j;

	//if (n < 0) nrerror("Negative factorial in routine factrl");
	if (n > 32)
	{
		// Larger value than size of table is required so use
		// the gamma function to evaluate it. 
		// Note it may overflow
		return exp(gammln(n+1.0));
	}

	while (ntop<n)
	{
		// Fill in table up to desired value.
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[n];
}


double plgndr(int l, int m, double x)
{
	// Computes the associated Legrendre polynomials
	// P(m,l){x}. Here m and l are integers satisfying
	// 0<=m<=l, while x lies in the range -1<=x<=1
	int i,ll;
	double fact,pll,pmm,pmmp1,somx2;

	double normalize = pow((2.0*l+1.0)/2.0*factrl(l-m)/factrl(l+m) ,0.5);

	// Compute P(m,m)
	pmm=1.0;
	if (m>0)
	{
		somx2=sqrt((1.0-x)*(1.0+x));
		fact=1.0;
		for (i=1;i<=m;i++)
		{
			pmm *= -fact*somx2;
			fact += 2.0;
		}
	}
	if (l==m)
	{
		return normalize*pmm;
	}
	else
	{
		// Compute P(m,m+1)
		pmmp1=x*(2*m+1)*pmm;
		if (l==(m+1))
		{
			return normalize*pmmp1;
		}
		else
		{
			// Compute P(m,l), l>m+1
			for (ll=m+2;ll<=l;ll++)
			{
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm=pmmp1;
				pmmp1=pll;
			}
			return normalize*pll;
		}
	}
}


double heval(double *x, double eta, double phi, double w, double a, double g, double Omega, double h0, int M, int kappa)
{
	double value = 0.0;
	int m,l,count;
	// The value of pi
	const double pi = acos(-1.0);
	double factor=sqrt(1.0/(2.0*pi));
	for (m=1;m<=M;m++)
	{
		for (l=kappa*m;l<=2*M+kappa*m-1;l++)
		{
			if (((l-m*kappa)%2)==0)
			{
				count = (m-1)*M+(l-kappa*m)/2;
				//value += x[2*M*M+(m-1)*M+(l-kappa*m)/2]*P[m][l-kappa*m][j]*C[m][i];
				value += *(x+count)*plgndr(l,kappa*m,sin(phi))*cos(kappa*m*eta)*factor;
			}
		}
	}

	value += h0+w*a*a*(2.0*Omega+w)/(2.0*g)*pow(cos(phi),2.0);

	return value;

}


void fsfield(double w, double a, double g, double Omega, int kappa, double *eta, double *phi, double *C, double *y, int mm, int nn, double h0, int M)
{
	// Computational function that calculates the pressure
	// at a series of grid points in latitude and longitude.
	int ii,jj;

	for (ii=0; ii<mm; ii++) {
		for (jj=0; jj<nn; jj++) 
		{
			//*(y+jj*mm+ii) = hz(w,a,g,Omega, *(phi+ii),h0)+epsilon*cos( *(eta+jj)*kappa)*h1(C,Clen,kappa, *(phi+ii));
			*(y+jj*mm+ii) = heval(C,*(eta+jj),*(phi+ii),w,a,g,Omega,h0,M,kappa);
		}
	}
  
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *eta, *phi, *C, *y;
  int M,i,j,kappa,etarows,etacols,phirows,phicols,mrows,ncols;
  double w,a,g,Omega,h0;
  
  /*  check for proper number of arguments */

  if(nrhs!=10) 
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
    mexErrMsgTxt("Input C must be a noncomplex double vector.");
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

  /* check to make sure the tenth input argument is a scalar */
  if( !mxIsDouble(prhs[9]) || mxIsComplex(prhs[9]) ||
      mxGetN(prhs[9])*mxGetM(prhs[9])!=1 ) {
    mexErrMsgTxt("Input M must be a scalar(integer).");
  }

  /*  get the scalar inputs Fr, M, Ro, Sr, Ka, w, kappa and epsilon, h0*/
  kappa = mxGetScalar(prhs[3]);
  w = mxGetScalar(prhs[4]);
  a = mxGetScalar(prhs[5]);
  g = mxGetScalar(prhs[6]);
  Omega = mxGetScalar(prhs[7]);
  h0 = mxGetScalar(prhs[8]);
  M = mxGetScalar(prhs[9]);

  /*  create a pointer to the input vectors eta, phi and C */
  eta = mxGetPr(prhs[0]);
  phi = mxGetPr(prhs[1]);
  C = mxGetPr(prhs[2]);

  /*  get the dimensions of the vectors eta, phi */
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
  fsfield(w,a,g,Omega,kappa,eta,phi,C,y,mrows,ncols,h0,M); 
}
