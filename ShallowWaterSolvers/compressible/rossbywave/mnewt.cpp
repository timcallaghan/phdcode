#include "nr.h"
#include "parameters.h"
//#include <stdlib.h>

using namespace std;

void NR::mnewt(Vec_IO_DP &coeffs,Vec_I_DP &phi,Mat_I_DP &Ceta,Mat_I_DP &Seta,Mat_I_DP &Cphi1,Mat_I_DP &Sphi1,Mat_I_DP &Cphi2,Mat_I_DP &Sphi2)
/* Given and initial guess x[0..n-1] for a root in n dimensions, take
ntrial Newton-Raphson steps to improve the root. Stop if the root
converges in either summed absolute variabe increments tolx or summed
absolute function values tolf. */
{
	int k,i,n,halvits;
	n = numunknowns;
	
	// Vector and Matrix objects used in the method
	DP errx,errf,errxtemp;
	Vec_DP p(n), fvec(n), delx(n), delxnew(n+1), fnew(n), xtemp(n), ftemp(n), tempdelx(n);
	Mat_DP fjac(n,n);

	// Used in QR factorisation
	bool singcheck;
	Vec_DP c(n);
	Vec_DP d(n);

	// Input and output file objects
	std::fstream fout;
	std::fstream rout;
	std::fstream lout;
	// Set output precision
	fout.precision(decimals);
	rout.precision(decimals);
	lout.precision(decimals);

	// Open a new file for output
	fout.open("coeffs.txt", std::ios::out);
	// Open a new file for the residual vector
	rout.open("residual.txt", std::ios::out);

	for (k=1;k<=NTRIAL;k++)
	{
		cout << "Starting iteration " << k << endl;

		if (k==1)
		{
			// write the initial values of the
			// coefficient vector to the file
			for (i=0;i<=n;i++)
			{
				if(i<n)
					fout << coeffs[i] << " ";
				else
					fout << coeffs[i] << std::endl;
			}
		}
		
		// Calculate the residual vector.
		residual(coeffs,phi,Ceta,Seta,Cphi1,Sphi1,Cphi2,Sphi2,fvec);
		cout << " Done the calculation of the residual" << endl;

		// Write the current value of the
		// residual vector to a file for
		// analysis.
		for (i=1;i<=n;i++)
		{
			if (i<n)
				rout << fvec[i-1] << " ";
			else
				rout << fvec[i-1] << std::endl;
		}

		// Check function convergence
		errf=0.0;
		for (i=0;i<n;i++) errf += fabs(fvec[i]);
		if (errf <= TOLF)
		{
			cout << "TOLF reached...exiting" << endl;
			fout.close();
			rout.close();
			return;
		}
		cout << "errf=" << errf << endl;
		
		// Evaluate the jacobian...since no convergence...
		jacobian(coeffs,phi,Ceta,Seta,Cphi1,Sphi1,Cphi2,Sphi2,fjac);
		cout << " Done the calculation of the jacobian" << endl;

		// Right hand side of linear equations.
		p = -1.0*fvec;

		// Solve the system using QR factorisation
		NR::qrdcmp(fjac,c,d,singcheck);
		NR::qrsolv(fjac,c,d,p);

	    // Translate this to a vector with all coeffs
		for (i=0;i<=n;i++)
		{
			if (i<fixedval)
			{
				delxnew[i]=p[i];
			}
			else
			{
				if (i>fixedval)
				{
					delxnew[i]=p[i-1];
				}
				else
				{
					// Insert a zero at the position of fixed
					delxnew[i]=0.0;
				}
			}
		}

		// Implement the halving technique
		// First calulate xtemp
		xtemp = coeffs + delxnew;
		// Now calculate the residual for this updated ftemp value
		residual(xtemp,phi,Ceta,Seta,Cphi1,Sphi1,Cphi2,Sphi2,fvec);
		// and the one-norm error
		errxtemp = 0.0;
		for (i=0;i<n;i++) errxtemp += fabs(fvec[i]);
		// Initialize halvings counter
		halvits=0;
		// Check to see if we need to use halving and if so do it
		while ((errxtemp > errf) && (halvits <= maxhalvings))
		{
			cout << "in here" << endl;
			delxnew = 0.5*delxnew;
			xtemp = coeffs + delxnew;
			residual(xtemp,phi,Ceta,Seta,Cphi1,Sphi1,Cphi2,Sphi2,fvec);
			errxtemp = 0.0;
			for (i=0;i<n;i++) errxtemp += fabs(fvec[i]);
			halvits++;
		}
		
		cout << "The number of halvings was: " << halvits << endl;
		// If maximum halving iterations was reached then we need to exit with an error
		if (halvits == maxhalvings+1)
		{
			NR::nrerror("Maximum halvings reached...start from new guess");
			return;
		}


		// Update the coeffs vector 
		coeffs = coeffs + delxnew;
		// Write the coefficients to file
		for (i=0;i<=n;i++)
		{
			if(i<n)
				fout << coeffs[i] << " ";
			else
				fout << coeffs[i] << std::endl;
		}
		// overwite the coefficients to another file
		// for use in another call to the program
		// using previous coefficient values.
		lout.open("current.txt", std::ios::out);
		// Write the coefficients to the file
		for (i=0;i<=n;i++)
		{
			if(i<n)
				lout << coeffs[i] << " ";
			else
				lout << coeffs[i] << std::endl;
		}
		lout.close();

		// Check for convergence
		errx = 0.0;
		for (i=0;i<=n;i++) errx += fabs(delxnew[i]);
		if (errx <= TOLX)
		{
			cout << "TOLX reached...exiting" << endl;
			fout.close();
			rout.close();
			return;
		}
		cout << "errx=" << errx << endl;
	}
	// Close the files
	fout.close();
	rout.close();
	return;
}

