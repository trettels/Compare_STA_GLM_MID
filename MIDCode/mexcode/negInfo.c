#include "mex.h"
#include <stdio.h>
#include <math.h>

/*
 *  negInf.c
 *  
 *
 *  Created by Sean Trettel on 6/28/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 *	Usage:
 *		[negInfo] = negInfo(px, pxt, bins)
 *
 *	Input:
 *		px: P(x) distribution
 *		pxt: P(x|spike) distribution
 *		bins: number of histogram bins
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs!=3) mexErrMsgTxt("Wrong number of input arguments"); //Error message thrown if there are an incorrect number of input arguments
	if(nlhs!=1) mexErrMsgTxt("Wrong number of output arguments"); //Error message thrown if there are an incorrect number of output arguments
	// Initialize varaibles
	double *px = mxGetPr(prhs[0]);
	double *pxt = mxGetPr(prhs[1]);
	int leng = (int)(mxGetScalar(prhs[2]));
	int i;
	double negInfo=0.0;

	for(i=0; i< leng; i++)
	{
		if(pxt[i]>0 && px[i] > 0)	//Only evaluate if P(x|spike) is nonzero
		{
			negInfo -=pxt[i]*log(pxt[i]/px[i]);  //Subtract this bin's contribution to the overall Kullback-Leibler divergence
		}
	}
	plhs[0]=mxCreateDoubleScalar(negInfo);  //Assign negInfo to the output structure
	
}
