#include "nr.h"
#include "fieldvars.h"
#include "dphysics.h"
#include "parameters.h"

using namespace std;


void NR::jacobian(Vec_IO_DP &coeffs,Vec_I_DP &phi,Mat_I_DP &Ceta,Mat_I_DP &Seta,Mat_I_DP &Cphi1,Mat_I_DP &Sphi1,Mat_I_DP &Cphi2, Mat_I_DP &Sphi2,Mat_IO_DP &jac)
/* Computes the jacobian matrix using
analytical expressions for the entries. This should
be correct to machine precision.
*/
{
	DP phival,wavespeed;
	Vec_DP vars(9);
	int i,j;

	wavespeed=coeffs[numunknowns];

	for (i=1;i<=M;i++)
	{
		// This loop corresponds to eta
		for (j=1;j<=N;j++)
		{
			// This loop corresponds to phi
			// Get the phi grid point
			phival = phi[j-1];
			// Set up the field variables vector
			vars[0] = ulameval(coeffs,i,j,phival,Ceta,Cphi1);
			vars[1] = uphieval(coeffs,i,j,Seta,Sphi2);
			vars[2] = heval(coeffs,i,j,Ceta,Cphi2);
			// eta derivatives
			vars[3] = ulamevaleta(coeffs,i,j,Seta,Cphi1);
			vars[4] = uphievaleta(coeffs,i,j,Ceta,Sphi2);
			vars[5] = hevaleta(coeffs,i,j,Seta,Cphi2);
			// phi derivatives
			vars[6] = ulamevalphi(coeffs,i,j,phival,Ceta,Sphi1);
			vars[7] = uphievalphi(coeffs,i,j,Seta,Cphi2);
			vars[8] = hevalphi(coeffs,i,j,Ceta,Sphi2);
    
			// Calculate the function derivatives
			// at the current colloction point in 
			// a row-wise manner.
			dmass(wavespeed,vars,phival,i,j,Ceta,Seta,Cphi1,Sphi1,Cphi2,Sphi2,jac);
			dlammomen(wavespeed,vars,phival,i,j,Ceta,Seta,Cphi1,Sphi1,Cphi2,Sphi2,jac);
			dphimomen(wavespeed,vars,phival,i,j,Ceta,Seta,Cphi1,Sphi1,Cphi2,Sphi2,jac);
		}
	}



	/*
	//////////////////////////////////////////////////////
	// Use this if we are going to split the integrals up
	//////////////////////////////////////////////////////
	// Calculate the Jacobian elements for the volume specification equation
	// Set up some memory
	Vec_DP yall(0.0,numunknowns+1);
	// Here we define the jacobian integral tolerance
	DP jacintol=1.0e-02;
	// Set up some temp variables to store results
	DP I1val, I2val, I3val;
	// First do the terms with m=0
	int mval=0;
	int nval=0;
	for (nval=0;nval<=N;nval++)
	{
		// These are the individual integral components when m=0.
		I1val = (2.0*gamma-1.0)/(gamma-1.0)*pow(gamma*Bval,2.0)*quad2dmod(I5,0.0,pi/kappa,0.0,pi/2.0,coeffs,jacintol,mval,nval);
		I2val = 2.0*(gamma-1.0)*gamma*Bval/arad*quad2dmod(I6,0.0,pi/kappa,0.0,pi/2.0,coeffs,jacintol,mval,nval);
		I3val = 2.0*pow((gamma-1.0)/arad,2.0)*quad2dmod(I7,0.0,pi/kappa,0.0,pi/2.0,coeffs,jacintol,mval,nval);
		yall[2*M*N+nval] = -1.0*Mbp*(I1val + I2val + I3val);
	}
	// Now take care of m>0 case
	for (mval=1;mval<=M-1;mval++)
	{
		for (nval=1;nval<=N;nval++)
		{
			I1val = (2.0*gamma-1.0)/(gamma-1.0)*pow(gamma*Bval,2.0)*quad2dmod(I5,0.0,pi/kappa,0.0,pi/2.0,coeffs,jacintol,mval,nval);
			I2val = 2.0*(gamma-1.0)*gamma*Bval/arad*quad2dmod(I6,0.0,pi/kappa,0.0,pi/2.0,coeffs,jacintol,mval,nval);
			I3val = 2.0*pow((gamma-1.0)/arad,2.0)*quad2dmod(I7,0.0,pi/kappa,0.0,pi/2.0,coeffs,jacintol,mval,nval);
			yall[2*M*N+mval*N+nval] = -1.0*Mbp*(I1val + I2val + I3val);
		}
	}
	*/



	//*
	//////////////////////////////////////////////////////
	// Use this if we don't want to split the integrals up
	//////////////////////////////////////////////////////
	// Calculate the Jacobian elements for the volume specification equation
	// Set up some memory
	Vec_DP yall(0.0,numunknowns+1);
	// Here we define the jacobian integral tolerance
	DP jacintol=1.0e-02;
	// First do the terms with m=0
	int mval=0;
	int nval=0;
	for (nval=0;nval<=N;nval++)
	{
		// These are the individual integral components when m=0.
		yall[2*M*N+nval] = -1.0*Mbp*quad2dmod(I4,0.0,pi/kappa,0.0,pi/2.0,coeffs,jacintol,mval,nval);
	}
	// Now take care of m>0 case
	for (mval=1;mval<=M-1;mval++)
	{
		for (nval=1;nval<=N;nval++)
		{
			yall[2*M*N+mval*N+nval] = -1.0*Mbp*quad2dmod(I4,0.0,pi/kappa,0.0,pi/2.0,coeffs,jacintol,mval,nval);
		}
	}
	//*/

    // Now fill out jac without translating the fixed coeff derivative
    for (int z=0;z<=numunknowns;z++)
	{
        if (z<fixedval)
		{
			jac[numunknowns-1][z]=yall[z];
		}
        else
		{
			if (z>fixedval)
			{
				jac[numunknowns-1][z-1]=yall[z];
			}
			else
			{
				// Do nothing since we don't want to translate the fixed coeff
			}
		}
	}
	// Now write the jacobian to a file for checking...
	// Input and output file objects
	//std::fstream fout;
	//fout.precision(decimals);
	// Open a new file for output
	//fout.open("jacobian.txt", std::ios::out);
	//fout << jac << endl;
	//fout.close();
	
	return;
}
