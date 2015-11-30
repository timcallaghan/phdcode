// File to store all the variable parameters
// that can be used in the problem.

// Prevent multiple inclusions of this
// header file...
#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <math.h>
#include "nr.h"
using namespace std;


// Declare static constants...
const DP pi = acos(-1.0); // The value of pi
const DP g = 9.80616; // Gravitational acceleration
const DP Omega = 2.0*pi/(24.0*60.0*60.0); // Earth's angular velocity
const DP aref = 6.37122e06; // Radius of the Earth
// Declare characteristic scales
const DP href = 8.0e03; // Characteristic length scale (m)
const DP cref_ = Omega/30.0; // Characteristic angular velocity scale of the Rossby wave (s^-1)
const DP vref = 40.0; // Characteristic velocity scale(m*s^-1)
const DP wref = 7.848e-06; // Characteristic angular velocity scale of the zonal flow (s^-1)
// Define the nondimensional numbers
const DP Sr = aref*cref_/vref; // Strouhal number
const DP Fr = vref/sqrt(g*href); // Froude number
const DP Ro = vref/(2.0*Omega*aref); // Rossby number
// Dimensionless zonal flow parameters
const DP arad = aref/href; // Dimensionless value of the Earth's radius...relative to href...

// Max output decimal precision
const int decimals=20;

// Tolerances and stopping criteria
const DP intol = 1.0e-14; // Quadrature integration tolerance
const DP TOLX = 1.0e-12; // Stopping criteria on delx
const DP TOLF = 1.0e-12; // Stopping criteria on f
const int maxhalvings = 30; // Upper limit on number of halvings allowed
const int NTRIAL = 30; // Maximum number of Newton steps to take....hopefully less than this!


// Define a class to load all the linearised constants.
class cLinParams
{
// Define the constants used in the calculation process.
public: 
	// The truncation level from the linear solution.
	int Nlin;
	// The number of wavelengths around a latitude circle
	int kappa;
	// The polar free surface height parameter
	DP h0;
	// The zonal angular velocity
	DP w;
	// The volume of the associate zonal flow for w and h0
	// Note that it's the hemispherical volume.
	DP Vzon;
	// Default constructor...
	cLinParams()
	{
		// Read in the dynamic values from file...
		ifstream inParams("linparams.txt", ios::in);
		if (!inParams)
		{
			cerr << "File linparams.txt could not be opened" << endl ;
			exit(1);
		}
		while (inParams >> kappa >> Nlin >> w >> h0 >> Vzon)
		inParams.close();
	}
};

// Define a class to load the old parameters.
class cOldParams
{
// Define the constants used in the calculation process.
public: 
	// Previous longitudinal truncation level.
	int Mold;
	// Previous latitudinal truncation level.
	int Nold;
	// Default constructor...
	cOldParams()
	{
		// Read in the dynamic values from file...
		ifstream inParams("oldparams.txt", ios::in);
		if (!inParams)
		{
			cerr << "File oldparams.txt could not be opened" << endl ;
			exit(1);
		}
		while (inParams >> Mold >> Nold)
		inParams.close();
	}
};

// Initialise one instance each of the classes
// and then assign the dynamic intialisation
// variables to their constant couterparts. 
const cLinParams data_;
const cOldParams odata;

// The number of coeffs in each series from the linear solution
const int Nlin = data_.Nlin;
// The number of wavelengths used in the linear solution.
const int kappa = data_.kappa;
// The polar free surface height used in the linear solution
const DP h0 = data_.h0;
// The zonal angular velocity used in the linear solution
const DP w = data_.w*wref*aref/vref;
// The RH amplitude condition used in the linera solution
const DP Vzon = data_.Vzon;

// Previous longitudinal truncation
const int Mold = odata.Mold;
// Previous latitudinal truncation
const int Nold = odata.Nold; 


// Define another class to hold all the dynamic constants.
class cNewParams
{
// Define the constants used in the calculation process.
public: 
	// Which type of bootstrapping we want to use
	int bootstrap;
	// The scaling factor for the amplitude
	DP scale;
	// New longitudinal truncation level.
	int M;
	// New latitudinal truncation level.
	int N;
	// The type of fixed coefficient
	int fixtype;
	// The index of the fixed coefficient
	int fixedval;
	// Default constructor...
	cNewParams()
	{
		// Read in the dynamic values from file...
		ifstream inParams("newparams.txt", ios::in);
		if (!inParams)
		{
			cerr << "File newparams.txt could not be opened" << endl ;
			exit(1);
		}
		while (inParams >> bootstrap >> scale >> fixtype)
		inParams.close();
		// Test for valid bootstrap values and assign M and N
		// appropriately
		switch(bootstrap)
		{
		case 1:
			M = Nlin;
			N = Nlin;
			break;
		case 0:
			M = Mold+1;
			N = Nold+1;
			break;
		case 2:
			M = Mold;
			N = Nold;
			break;
		default:
			cerr << "Illegal bootstrap value" << endl ;
			exit(1);
		}
		// Check to see if we have a valid type of fixed coeff
		switch(fixtype)
		{
		case 0:
			// Holding the first free surface coeff fixed (not in zonal flow part)
			fixedval = 2*M*N+N+1;
			break;
		case 1:
			// Holding the wavespeed fixed...
			fixedval = 3*M*N+1;
			break;
		default:
			cerr << "Illegal value for fixtype in newparams.txt" << endl ;
			exit(1);
		}
	}
};


const cNewParams ndata;
// Value of bootstrap
const int bootstrap = ndata.bootstrap;
const int M = ndata.M;
const int N = ndata.N;
// Index of fixed coefficient...first free surface coeff (not in zonal flow part)
const int fixedval = ndata.fixedval; 
// The number of unknowns in the nonlinear model
const int numunknowns = 3*M*N+1;
// The amplitude scaling factor
const DP scale = ndata.scale;

#endif
