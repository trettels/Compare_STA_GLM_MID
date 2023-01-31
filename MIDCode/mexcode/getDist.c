#include "mex.h"
#include <stdio.h>
#include <math.h>
/*
 *  getDist.h
 *  
 *
 *  Created by Sean Trettel on 6/28/10.
 *  Fairhall Lab, Department of Physiology and Biophysics, University of Washington.
 *	
 *	This function creates the P(x) and P(x|spike) distributions, given the spike train and projected stimulus
 *	
 *	USAGE:
 *		[px, pxt] = getDist(nSpikes, nTrials, min, step, BINS, nfilt, stimProj, spikes);
 *
 *	INPUTS:
 *		prhs[0] = Number of Spikes
 *		prhs[1] = trial length (same as length(stimProj))
 *		prhs[2] = minimum projection value (center of the leftmost bin)
 *		prhs[3] = Step size (distance between histogram bin centers)
 *		prhs[4] = Number of bins
 *		prhs[5] = length of the filter, used to compensate for spike offset
 *		prhs[6] = Stimulus projection vector
 *		prhs[7] = Spike train
 *
 *	OUTPUTS:
 *		plhs[0] = P(x)
 *		plhs[1] = P(x|spike)
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs!=9) mexErrMsgTxt("Wrong number of input arguments"); //Error message thrown if there are an incorrect number of input arguments
	if(nlhs!=2) mexErrMsgTxt("Wrong number of output arguments"); //Error message thrown if there are an incorrect number of output arguments
	
	// Set up the variables
	int index, cross_index;
	int j,k;
	int nSpikes = (int)(mxGetScalar(prhs[0]));
	int nTrials = (int)(mxGetScalar(prhs[1]));
	double r=1/(double)(nTrials);
	//double min= mxGetScalar(prhs[2]);
	//double step=mxGetScalar(prhs[3]);
    mxArray *xMin=mxDuplicateArray(prhs[2]);
    mxArray *xStep = mxDuplicateArray(prhs[3]);
	int BINS=(int)(mxGetScalar(prhs[4]));
	int nfilt=(int)(mxGetScalar(prhs[5]));
	mxArray *stimProj=mxDuplicateArray(prhs[6]);
	double *spikes = mxGetPr(prhs[7]);
	int nvec = (int)(mxGetScalar(prhs[8]));
	
	// Set up the output structures
    int distsize = (int)(pow(BINS,nvec));
	plhs[0]=mxCreateDoubleMatrix(1,distsize,mxREAL);
	plhs[1]=mxCreateDoubleMatrix(1,distsize,mxREAL);
	double *px=mxGetPr(plhs[0]);
	double *pxt=mxGetPr(plhs[1]);
//	mexPrintf("\nNumSpikes: %u\n",nSpikes);
/*
    for(k=0;k<nTrials; k++)
    {
        if(spikes[k]>0) mexPrintf("\nIndex: %u, spike:%f\n",k,spikes[k]);
    }
*/
	for(j=0; j<nTrials-1; j++)
	{
        int cross_index[nvec];
        for(k=0;k<nvec;k++) cross_index[k]=0;
	    //double *stp = mxGetPr(stimProj);
        //double min = mxGetScalar(xMin);
	    //double step = mxGetScalar(xStep);
        for(k=0; k<nvec; k++)
        {
            // load current directions' values
            //mxArray *st = mxGetCell(stimProj, k);
            double *stp=mxGetPr(mxGetCell(stimProj, k));
            double min=mxGetScalar(mxGetCell(xMin,k));
            double step=mxGetScalar(mxGetCell(xStep,k));
            // Assign the cross-index
            cross_index[k]=(int)(round((stp[j]-min)/step));
            //mexPrintf("\nBins: %u,  index: %7.7f,  stp:  %7.7f,  step: %f  \n",BINS, cross_index[k],stp[j], step);
            //mexPrintf("\nmin:  %f, step: %f, index: %f\n",min,step,cross_index[k]);
            if(cross_index[k] >= BINS)
            {
				cross_index[k]=BINS-1;  
				//mexPrintf("+");
			}
            else if(cross_index[k] < 0)
            {
                	cross_index[k]=0;
                //mexPrintf("-");
			}
            //if(cross_index[k]==BINS) cross_index[k]=BINS-1;
            //else if(cross_index[k] < 0){ mexPrintf("\nindex: %f",cross_index[k]); mexPrintf("\nk: %i",k);  mexErrMsgTxt("cross_index Index Out Of Bounds (too small) in getDist.c"); }
            //else if(cross_index[k] >BINS){ mexPrintf("\nindex: %f",cross_index[k]);  mexPrintf("\nk: %i",k);  mexErrMsgTxt("cross_index Index Out Of Bounds (too large) in getDist.c"); }
        }
        // Assign the complete index for the probability array
        index=0;
        for(k=0;k<nvec;k++) index += (int)(round(cross_index[k]*pow(BINS,k)));  //*(pow(BINS,k))));

		px[index] += 1/(double)(nTrials);					//Bin the P(x) value
		pxt[index]+= spikes[j+nfilt]/(double)(nSpikes);		//Bin the P(x|spike) value, if there is a spike present at the end of the window  //spikes[j+nfilt]/(double)(nSpikes);
		//mexPrintf("\nIndex: %u, relevant spike: %f",index,spikes[j+nfilt]);
	}
	
	//Free Memory
	
	
}
