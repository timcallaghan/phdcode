/*=============================================================
 * makealpha.c 
 * Provides a quick way of generating an alpha channel for
 * an PNG file. The result is double matrix with
 * alpha channel information (1.0=clear, 0.0=opaque) 
 *
 * This is a MEX-file for MATLAB.  
 *============================================================*/

/* $ Revision: 1.5 $ */
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* Declare variables. */ 
    int picsize, i;
	unsigned char *pRed, *pGreen, *pBlue;
	unsigned char toll;
    double *y;
	double clearval, opaqueval;
    
    /* Check for proper number of input and output arguments. */    
    if (nrhs != 2) {
    mexErrMsgTxt("Two input arguments required.");
    } 
    if (nlhs > 1){
    mexErrMsgTxt("Too many output arguments.");
    }
    
	// Get the tollerance...need to cast it to an unsigned char
	toll = (unsigned char) mxGetScalar(prhs[1]);
    
    // Get the pointer to the first data entry in the array.
	// We define this to be the Red colour channel. Note that the
	// PNG data is stored as uint8 in matlab (so unsigned char in c)
    pRed = (unsigned char *)mxGetPr(prhs[0]);

	// Get the size of the array (it's square so only need 1 dimension)
	picsize = mxGetM(prhs[0]);

	// Now set the pointers to the Green and Blue channels
	// We just point to a different part of prhs[0]...note that
	// the memory addition works relative to the type of the pointer
	// ...in this case an unsigned char...
	pGreen = pRed + picsize*picsize;
	pBlue = pRed + 2*picsize*picsize;

	/*  set the output pointer to the output matrix */
	plhs[0] = mxCreateDoubleMatrix(picsize,picsize,mxREAL);
  
	/*  create a C pointer to a copy of the output matrix */
	y = mxGetPr(plhs[0]);

	// Now we should have pointers to the three colour channels
	// and also the output matrix so we no longer need to treat 
	// the pointer to the RHS as multidimensional...just treat
	// like normal matrices

	clearval = 1.0;
	opaqueval = 0.0;

	// Just loop through and check every element to see if the tollerance
	// is satisfied...
	for (i=0;i<picsize*picsize;i++)
	{
		// First check if the red component is greater than or 
		// equal to the tolerance
		if ( *(pRed+i) >= toll)
		{
			// We've passed the first test...now check the green channel
			if ( *(pGreen+i) >= toll)
			{
				// We've passed the second test...now test the blue channel
				if ( *(pBlue+i) >= toll)
				{
					// All conditions are met so we make the pixel transparent
					*(y+i) = clearval;
				}
				else
				{
					// We have a darked colour so set it to be opaque
					*(y+i) = opaqueval;
				}
			}
			else
			{
				// We have a darked colour so set it to be opaque
				*(y+i) = opaqueval;
			}
		}
		else
		{
			// We have a darked colour so set it to be opaque
			*(y+i) = opaqueval;
		}
	}
}

