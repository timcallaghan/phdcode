// This program attempts to solve the full nonlinear
// shallow atmosphere equations on the suface of a 
// sphere. The field variables are approximated by 
// truncated Fourier series and the governing 
// equations for the flow field are then 
// evaluated at a set of discrete internal nodes, 
// using the method of collocation. This whole process
// is then iterated using Newton's method until 
// convergence occurs and the residual is sufficiently
// small enough.
//
// Some points to note...
//
// 1) The Fourier series approximations are constructed so
//    that the fluid velocity is zero at the poles. In addition
//    it is assumed that the latitudinal velocity component is
//    anit-symmetric with respect to the equator whilst all
//    other field variables are assumed symmetric with respect
//    to the equator.
// 2) The number of wavelengths around a latitude circle can
//    be specified in advance. This way the computational
//    domain can be reduced in size for greater accuracy
//    over the whole sphere.
// 3) The governing equations involve time dependence but
//    only in the combination eta=lambda-ct where c is the
//    speed of the travelling progresive wave. 	  

#include "nr.h"
#include "parameters.h"
#include "physics.h"
#include "fieldvars.h"

using namespace std;

int main(void)
{
	// Declare a counting variable
	int i;

	// Set the output precision
	cout.precision(decimals);

	// Calculate the collocation points in phi
	DP delphi = pi/(2.0*N);
	DP epsphi = delphi/2.0;
	Vec_DP phi(N);
	for (i=1;i<=N;i++)
	{
		phi[i-1]=(i-1)*delphi+epsphi;
	}
	// Now cache the basis functions for efficient look-up.
	Mat_DP Ceta(M,M),Seta(M,M),Cphi1(N+1,N),Sphi1(N+1,N),Cphi2(N+1,N),Sphi2(N+1,N);
	NR::cacheCeta(Ceta);
	NR::cacheSeta(Seta);
	NR::cacheCphi1(Cphi1,phi);
	NR::cacheSphi1(Sphi1,phi);
	NR::cacheCphi2(Cphi2,phi);
	NR::cacheSphi2(Sphi2,phi);

	// Initialise a starting guess for the
	// vector of unknown coefficients to zero...note we also include
	// the fixed coeff in this vector
	Vec_DP coeffs(0.0,numunknowns+1);

	// We now read in the starting guess for the 
	// coefficients from file using our routines.
	// Implement bootstrap method based on switch function
	switch(bootstrap)
	{
	case 0:
		// Load the previous guess with M and N changing
		// but amplitude is constant
		NR::loadNLinCoeffs(coeffs);
		break;
	case 1:
		// Load the linearised coeffs
		NR::loadLinCoeffs(coeffs);
		break;
	case 2:
		// Load the previous guess with M and N constant
		// and amplitude changing by our scale factor
		NR::loadNLinCoeffs(coeffs);
		coeffs[fixedval] = scale*coeffs[fixedval];
		break;
	default:
		NR::nrerror("Bootstrap value is incorrect...must be 0,1 or 2...");
	}
	
	// We can now do the Newton iterations on the
	// intial guess by calling mnewt.
	NR::mnewt(coeffs,phi,Ceta,Seta,Cphi1,Sphi1,Cphi2,Sphi2);

	// Now we update the value of Mold and Nold in oldparams.txt
	// Input and output file objects
	fstream fout;
	// Set output precision
	fout.precision(decimals);
	// Open a new file for output
	fout.open("oldparams.txt", std::ios::out);
	fout << M << endl;
	fout << N << endl;
	fout.close();

	return 0;
}