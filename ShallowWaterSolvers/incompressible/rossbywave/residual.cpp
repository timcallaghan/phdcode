#include "nr.h"
#include "fieldvars.h"
#include "physics.h"
#include "parameters.h"

using namespace std;


/////////////////////
// The residual function that returns the residual vector
// for the equations of motion.
void NR::residual(Vec_IO_DP &coeffs,Vec_I_DP &phi,Mat_I_DP &Ceta,Mat_I_DP &Seta,Mat_I_DP &Cphi1,Mat_I_DP &Sphi1,Mat_I_DP &Cphi2, Mat_I_DP &Sphi2,Vec_IO_DP &fvec)
{
	DP phival,wavespeed,Vnl;
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
        
			// Calculate the function values
			// at the current colloction point.
			fvec[(i-1)*3*M+(j-1)*3] = mass(vars,phival,wavespeed);
			fvec[(i-1)*3*M+(j-1)*3+1] = lammomen(vars,phival,wavespeed);
			fvec[(i-1)*3*M+(j-1)*3+2] = phimomen(vars,phival,wavespeed);
		}
	}
	// Calculate the current volume for the free surface
	Vnl=2.0*kappa*arad*arad*quad2d(surfint,0.0,pi/kappa,0.0,pi/2.0,coeffs,intol);
	// Now put in the mass specification condition
	fvec[numunknowns-1]=1.0-Vnl/Vzon;
	return;
}
