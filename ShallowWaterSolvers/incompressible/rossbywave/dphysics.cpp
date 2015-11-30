#include "nr.h"
#include "dphysics.h"
#include "parameters.h"

// Functions that evaluate the derivatives
// of the equations of motion with respect to a coefficient
// at a given point.

void dmass(const DP wavespeed, Vec_I_DP &vars,const DP phival,const int i,const int j,Mat_I_DP &Ceta,Mat_I_DP &Seta,Mat_I_DP &Cphi1,Mat_I_DP &Sphi1,Mat_I_DP &Cphi2, Mat_I_DP & Sphi2,Mat_IO_DP &jac)
{
    // Initialise the vector that will store all the derivatives.
	// We calculate them all and then pick of the ones we want to
	// include...a bit wasteful but easier to code...
    Vec_DP yall(numunknowns+1);
    //Terms from using the chain rule on the mass equation
    DP A1 = vars[5];
	DP A2 = vars[2];
	DP B1 = cos(phival)*vars[8]-vars[2]*sin(phival);
	DP B3 = vars[2]*cos(phival);
	DP C1 = vars[3]+cos(phival)*vars[7]-sin(phival)*vars[1];
	DP C2 = vars[0]-Sr*wavespeed*cos(phival);
	DP C3 = cos(phival)*vars[1];
	DP D = -Sr*cos(phival)*vars[5];

    // First we need to calculate the dervative for the uncommon indices in each series.
    // That is, the indices that are not part of a double series expression.
    // Terms related to m=0
    // Set up derivative with respect to the constant term H00
    yall[2*M*N] = C1;
    // Now do the rest of the eta-independent terms
	int m,n;
    for (n=1;n<=N;n++)
	{
        yall[2*M*N+n] = C1*Cphi2[n-1][j-1]-C3*2.0*n*Sphi2[n-1][j-1];
	}
    // Terms related to m=M
    for (n=1;n<=N;n++)
	{
        yall[M*N-N+n-1] = A1*Ceta[M-1][i-1]*Cphi1[n-1][j-1]-A2*kappa*M*Seta[M-1][i-1]*Cphi1[n-1][j-1];
        yall[2*M*N-N+n-1] = B1*Seta[M-1][i-1]*Sphi2[n-1][j-1]+B3*2.0*n*Seta[M-1][i-1]*Cphi2[n-1][j-1];
	}

    // Now do the common series components
    for (m=1;m<=M-1;m++)
	{
        for (n=1;n<=N;n++)
		{
            if (n==1)
			{
                // We need to use an alternative method for the H(i,1) derivatives
                yall[(m-1)*N+n-1] = A1*Ceta[m-1][i-1]*Cphi1[n-1][j-1]-A2*kappa*m*Seta[m-1][i-1]*Cphi1[n-1][j-1];
                yall[M*N+(m-1)*N+n-1] = B1*Seta[m-1][i-1]*Sphi2[n-1][j-1]+B3*2.0*n*Seta[m-1][i-1]*Cphi2[n-1][j-1];
                yall[2*M*N+m*N+n] = -C1*Ceta[m-1][i-1]*(Cphi2[0][j-1]+1.0)+C2*kappa*m*Seta[m-1][i-1]*(Cphi2[0][j-1]+1.0)+C3*Ceta[m-1][i-1]*2.0*Sphi2[0][j-1];
			}
            else
			{
                yall[(m-1)*N+n-1] = A1*Ceta[m-1][i-1]*Cphi1[n-1][j-1]-A2*kappa*m*Seta[m-1][i-1]*Cphi1[n-1][j-1];
                yall[M*N+(m-1)*N+n-1] = B1*Seta[m-1][i-1]*Sphi2[n-1][j-1]+B3*2.0*n*Seta[m-1][i-1]*Cphi2[n-1][j-1];
                yall[2*M*N+m*N+n] = C1*Ceta[m-1][i-1]*pow(-1.0,(double) n)*(Cphi2[n-1][j-1]+Cphi2[n-2][j-1])-C2*kappa*m*Seta[m-1][i-1]*pow(-1.0,(double) n)*(Cphi2[n-1][j-1]+Cphi2[n-2][j-1])-C3*Ceta[m-1][i-1]*pow(-1.0,(double) n)*(2.0*n*Sphi2[n-1][j-1]+2.0*(n-1)*Sphi2[n-2][j-1]);
			}
		}
	}
    // Fill out the derivative with respect to the wavespeed c
	yall[numunknowns] = D;
    // Now fill out jac without translating the fixed coeff derivative
    for (int z=0;z<=numunknowns;z++)
	{
        if (z<fixedval)
		{
			jac[(i-1)*3*M+(j-1)*3][z]=yall[z];
		}
        else
		{
			if (z>fixedval)
			{
				jac[(i-1)*3*M+(j-1)*3][z-1]=yall[z];
			}
			else
			{
				// Do nothing since we don't want to translate the fixed coeff
			}
		}
	}
	return;
}


void dlammomen(const DP wavespeed, Vec_I_DP &vars,const DP phival,const int i,const int j,Mat_I_DP &Ceta,Mat_I_DP Seta,Mat_I_DP &Cphi1,Mat_I_DP &Sphi1,Mat_I_DP &Cphi2, Mat_I_DP & Sphi2,Mat_IO_DP &jac)
{
    // Initialise the vector that will store all the derivatives.
	// We calculate them all and then pick of the ones we want to
	// include...a bit wasteful but easier to code...
    Vec_DP yall(numunknowns+1);
    // Terms from using the chain rule on the lambda momentum equation
	DP A1 = vars[3]-sin(phival)*vars[1];
	DP A2 = vars[0]-Sr*wavespeed*cos(phival);
	DP A3 = cos(phival)*vars[1];
	DP B1 = cos(phival)*vars[6]-(cos(phival)/Ro+vars[0])*sin(phival);
	DP C2 = 1.0/(Fr*Fr);
	DP D = -Sr*cos(phival)*vars[3];
    
	// First we need to calculate the dervative for the uncommon indices in each series.
    // That is, the indices that are not part of a double series expression.
    // Terms related to m=0
    // Set up derivative with respect to the constant term H00
    yall[2*M*N]=0.0;
	int m,n;
    // Now do the rest of the eta-independent terms
    for (n=1;n<=N;n++)
	{
        yall[2*M*N+n]=0.0;
	}
    // Terms related to m=M
    for (n=1;n<=N;n++)
	{
        yall[M*N-N+n-1] = A1*Ceta[M-1][i-1]*Cphi1[n-1][j-1]-A2*kappa*M*Seta[M-1][i-1]*Cphi1[n-1][j-1]-A3*(2.0*n-1.0)*Ceta[M-1][i-1]*Sphi1[n-1][j-1];
        yall[2*M*N-N+n-1] = B1*Seta[M-1][i-1]*Sphi2[n-1][j-1];
	}
    // Now do the common series components
    for (m=1;m<=M-1;m++)
	{
        for (n=1;n<=N;n++)
		{
            if (n==1)
			{
                // We need to use an alternative method for the H(i,1) derivatives
                yall[(m-1)*N+n-1] = A1*Ceta[m-1][i-1]*Cphi1[n-1][j-1]-A2*kappa*m*Seta[m-1][i-1]*Cphi1[n-1][j-1]-A3*(2.0*n-1.0)*Ceta[m-1][i-1]*Sphi1[n-1][j-1];
                yall[M*N+(m-1)*N+n-1] = B1*Seta[m-1][i-1]*Sphi2[n-1][j-1];
                yall[2*M*N+m*N+n] = C2*kappa*m*Seta[m-1][i-1]*(Cphi2[0][j-1]+1.0);
            }
			else
			{
                yall[(m-1)*N+n-1] = A1*Ceta[m-1][i-1]*Cphi1[n-1][j-1]-A2*kappa*m*Seta[m-1][i-1]*Cphi1[n-1][j-1]-A3*(2.0*n-1.0)*Ceta[m-1][i-1]*Sphi1[n-1][j-1];
                yall[M*N+(m-1)*N+n-1] = B1*Seta[m-1][i-1]*Sphi2[n-1][j-1];
                yall[2*M*N+m*N+n] = -C2*kappa*m*Seta[m-1][i-1]*pow(-1.0,(double) n)*(Cphi2[n-1][j-1]+Cphi2[n-2][j-1]);
            }
		}
	}
    // Fill out the derivative with respect to the wavespeed c
	yall[numunknowns] = D;
    // Now fill out jac without translating the fixed coeff derivative
    for (int z=0;z<=numunknowns;z++)
	{
        if (z<fixedval)
		{
			jac[(i-1)*3*M+(j-1)*3+1][z]=yall[z];
		}
        else
		{
			if (z>fixedval)
			{
				jac[(i-1)*3*M+(j-1)*3+1][z-1]=yall[z];
			}
			else
			{
				// Do nothing since we don't want to translate the fixed coeff
			}
		}
	}
	return;
}


void dphimomen(const DP wavespeed, Vec_I_DP &vars,const DP phival,const int i,const int j,Mat_I_DP &Ceta,Mat_I_DP &Seta,Mat_I_DP &Cphi1,Mat_I_DP &Sphi1,Mat_I_DP &Cphi2, Mat_I_DP & Sphi2,Mat_IO_DP &jac)
{
    // Initialise the vector that will store all the derivatives.
	// We calculate them all and then pick off the ones we want to
	// include...a bit wasteful but easier to code...
    Vec_DP yall(numunknowns+1);
    // Terms from using the chain rule on the phi momentum equation
    DP A1 = vars[4]+(cos(phival)/Ro+2.0*vars[0])*sin(phival);
	DP B1 = cos(phival)*vars[7];
	DP B2 = vars[0]-Sr*wavespeed*cos(phival);
	DP B3 = cos(phival)*vars[1];
	DP C3 = cos(phival)/(Fr*Fr);
	DP D = -Sr*cos(phival)*vars[4];

    // First we need to calculate the dervative for the uncommon indices in each series.
    // That is, the indices that are not part of a double series expression.
    // Terms related to m=0
    // Set up derivative with respect to the constant term H00
    yall[2*M*N] = 0.0;
	int m,n;
    // Now do the rest of the eta-independent terms
    for (n=1;n<=N;n++)
	{
        yall[2*M*N+n] = -C3*2.0*n*Sphi2[n-1][j-1];
	}
    // Terms related to m=M
    for (n=1;n<=N;n++)
	{
        yall[M*N-N+n-1] = A1*Ceta[M-1][i-1]*Cphi1[n-1][j-1];
        yall[2*M*N-N+n-1] = B1*Seta[M-1][i-1]*Sphi2[n-1][j-1]+B2*kappa*M*Ceta[M-1][i-1]*Sphi2[n-1][j-1]+B3*2.0*n*Seta[M-1][i-1]*Cphi2[n-1][j-1];
    }
    // Now do the common series components
    for (m=1;m<=M-1;m++)
	{
        for (n=1;n<=N;n++)
		{
            if (n==1)
			{
                // We need to use an alternative method for the H(i,1) derivatives
                yall[(m-1)*N+n-1] = A1*Ceta[m-1][i-1]*Cphi1[n-1][j-1];
                yall[M*N+(m-1)*N+n-1] = B1*Seta[m-1][i-1]*Sphi2[n-1][j-1]+B2*kappa*m*Ceta[m-1][i-1]*Sphi2[n-1][j-1]+B3*2.0*n*Seta[m-1][i-1]*Cphi2[n-1][j-1];
                yall[2*M*N+m*N+n] = C3*Ceta[m-1][i-1]*2.0*Sphi2[0][j-1];
            }
			else
			{
                yall[(m-1)*N+n-1] = A1*Ceta[m-1][i-1]*Cphi1[n-1][j-1];
                yall[M*N+(m-1)*N+n-1] = B1*Seta[m-1][i-1]*Sphi2[n-1][j-1]+B2*kappa*m*Ceta[m-1][i-1]*Sphi2[n-1][j-1]+B3*2.0*n*Seta[m-1][i-1]*Cphi2[n-1][j-1];
                yall[2*M*N+m*N+n] = -C3*Ceta[m-1][i-1]*pow(-1.0,(double) n)*(2.0*n*Sphi2[n-1][j-1]+2.0*(n-1.0)*Sphi2[n-2][j-1]);
            }
		}
	}
    // Fill out the derivative with respect to the wavespeed c
	yall[numunknowns] = D;
    // Now fill out jac without translating the fixed coeff derivative
    for (int z=0;z<=numunknowns;z++)
	{
        if (z<fixedval)
		{
			jac[(i-1)*3*M+(j-1)*3+2][z]=yall[z];
		}
        else
		{
			if (z>fixedval)
			{
				jac[(i-1)*3*M+(j-1)*3+2][z-1]=yall[z];
			}
			else
			{
				// Do nothing since we don't want to translate the fixed coeff
			}
		}
	}
	return;
}
