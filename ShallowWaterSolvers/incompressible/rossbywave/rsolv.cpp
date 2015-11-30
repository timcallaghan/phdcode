#include "nr.h"

void NR::rsolv(Mat_I_DP &a, Vec_I_DP &d, Vec_IO_DP &b)
{
	// Solves the set of n linear equations R.x=b, where R is an
	// upper triangular matrix stored in a and d. a[0..n-1][0..n-1]
	// and d[0..n-1] are input as the output of the routine qrdcmp
	// and are not modified. b[0..n-1] is input as the right-hand
	// side vector, and is overwritten with the solution vector
	// on output.

	int i,j;
	DP sum;

	int n=a.nrows();
	b[n-1] /= d[n-1];
	for (i=n-2;i>=0;i--)
	{
		for (sum=0.0,j=i+1;j<n;j++)
		{
			sum += a[i][j]*b[j];
		}
		b[i] = (b[i]-sum)/d[i];
	}
}

