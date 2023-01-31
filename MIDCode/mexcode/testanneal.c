#include "mex.h"
#include <stdio.h>
#include <math.h>


/*
 *  Anneal.c
 *  
 *
 *  Created by Sean Trettel on 6/28/10.
 *  Fairhall Lab, Department of Physiology and Biophysics, University of Washington
 *
 *
 *	USAGE:
 *		[bestVec, bestEnergy] = anneal(coolCoef, initTemp, maxRej, maxSuccess, maxTry, maxIter, stopTemp, heatCoeff, minEnergy, testVec, stim, spikes, nfilt, bins)
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Errors
	if(nrhs!=22) mexErrMsgTxt("Wrong number of input arguments"); //Error message thrown if there are an incorrect number of input arguments
	if(nlhs!=2) mexErrMsgTxt("Wrong number of output arguments"); //Error message thrown if there are an incorrect number of output arguments
	
	// Initialize and assign input varaibles
	double coolCoef=mxGetScalar(prhs[0]);
	double initTemp = mxGetScalar(prhs[1]);
	int maxRej = (int)(mxGetScalar(prhs[2]));
	int maxSuccess = (int)(mxGetScalar(prhs[3]));
	int maxTry = (int)(mxGetScalar(prhs[4]));
	int maxIter = (int)(mxGetScalar(prhs[5]));
	double stopTemp = mxGetScalar(prhs[6]);
	double heatCoef = mxGetScalar(prhs[7]);
	double minEnergy = mxGetScalar(prhs[8]);
	
	mxArray *testV = mxDuplicateArray(prhs[9]);
	mxArray *stim = mxDuplicateArray(prhs[10]);
	mxArray *test_stim = mxDuplicateArray(prhs[11]);
	mxArray *spikes = mxDuplicateArray(prhs[12]);
	mxArray *test_spikes = mxDuplicateArray(prhs[13]);
	mxArray *nfilt = mxDuplicateArray(prhs[14]);
	mxArray *bins = mxDuplicateArray(prhs[15]);
	double bestEnergy = mxGetScalar(prhs[16]);
	double oldEnergy = mxGetScalar(prhs[17]);
	double maxTemp = mxGetScalar(prhs[18]);
	//const mwSize *dimensions=mxGetDimensions(prhs[9]);
	//mwSize ndim=mxGetNumberOfDimensions(prhs[9]);
	int nvec =(int)(mxGetScalar(prhs[21]));
	//double *bestVec = mxGetPr(prhs[19]);		prhsT2[4]=bins;
    //int k;
    //for(k=0;k<ndim;k++) dimensions[k]=(int)dimensions[k];

    //int HCURRENT = (int)(mxGetScalar(mexGetVariable("global","HCURRENT")));
    
	
	int length = (int)(mxGetScalar(prhs[20]));
	// Initialize and assign output variables
	//plhs[0]=mxDuplicateArray(prhs[19]);
	
	//mxArray *outVec = mxGetPr(plhs[0]);
	double *outEnergy;
	
	double newEnergy, testEnergy=0.0;
	
	int i,j,k;
	//int powell_tol;

	int iTry=0,success=0,consec=0;
	double T = initTemp;
	
	mxArray *bestV = mxDuplicateArray(prhs[19]);
	mxArray *currentV=mxDuplicateArray(prhs[9]);
	mxArray *parentV= mxDuplicateArray(prhs[9]);
	mxArray *newV = mxDuplicateArray(prhs[9]);
	plhs[0] = mxDuplicateArray(bestV);//outVec[j]=bestVec[j];
	plhs[1] = mxCreateDoubleScalar(T);
    /*
	double *bestVec = mxGetPr(bestV);
	double *currentVec = mxGetPr(currentV);
	double *parentVec = mxGetPr(parentV);
	double *newVec = mxGetPr(newV);
	
	mxArray *plhsT2[1], *prhsT2[7];
	mxArray *plhsT1[1], *prhsT1[6];
	mxArray *plhsT3[1], *prhsT3[8];
	mxArray *plhsT4[1], *prhsT4[6];
	mxArray *plhsT5[1], *prhsT5[8];
	mxArray *plhsT6[1], *prhsT6[6];
	mxArray *plhsT7[1], *prhsT7[0];
	mxArray *plhsT8[1], *prhsT8[6];
    */

    
	/*
	//Initialize all vectors to the initial guess.
	for(j=0; j<nvec; j++)
	{
		mxSetCell(currentV, j, mxGetCell(testV,j));//currentVec[j]=testVec[j];
		mxSetCell(parentV, j, mxGetCell(testV,j));//parentVec[j]=testVec[j];
		mxSetCell(newV, j, mxGetCell(testV,j));//newVec[j]=testVec[j];
	}
	*/	
	mexPrintf("\t\t Beginning Annealing Algorithm\n");
	for(i=0; i<maxIter; i++)
	{
		//for(j=0; j<nvec; j++) mxSetCell(currentV, j, mxGetCell(parentV,j));//currentVec[j]=parentVec[j];
        currentV=mxDuplicateArray(parentV);
		iTry++;
		
		if(consec >= maxRej)
		{
			heatCoef *= .95;
			currentV=mxDuplicateArray(bestV);
		}
		
		if(iTry >= maxTry || success >= maxSuccess)
		{
			T *= coolCoef;
			if(T< stopTemp) T=stopTemp;
			plhs[1] = mxCreateDoubleScalar(T);
			iTry=1;
			success=1;
		}
		
/*
 * Generator Function
 */
		// Generate and assign the input variable structures...

		/*
        for(k=0;k<nvec;k++)
        {
            double norm = 0.0;
            mxArray *cv=mxDuplicateArray(mxGetCell(currentV, k));
            mxArray *bv=mxDuplicateArray(mxGetCell(bestV,k));
            double *currentVec = mxGetPr(cv);
            double *bestVec = mxGetPr(bv);
            for(j=0;j<length;j++) norm += currentVec[j]*currentVec[j];
            norm = sqrt(norm);
            //mexPrintf("%7.5f \n",norm);
            if(norm > 10000000) {
                mexPrintf("renormalizing\n");
                double bnorm = 0.0;
                for(j=0;j<length;j++) bnorm += bestVec[j]*bestVec[j];
                bnorm = sqrt(bnorm);
                for(j=0;j<length;j++) currentVec[j] = bestVec[j] / bnorm;
                mexPrintf("\n Re-scaled\n ");
                mxSetCell(currentV, k, cv);
                mxSetCell(bestV, k, bv);
            }
        }
        */
        //powell_tol =(int)(floor(i/10));
        	mxArray *prhsT1[7];
        	mxArray *plhsT1[1];
		int nlhsT=1, nrhsT=7;
		prhsT1[0]=currentV;
		prhsT1[1]=stim; //stim
		prhsT1[2]=spikes; //spikes
		prhsT1[3]=nfilt; //nfilt
		prhsT1[4]=bins; //bins
		prhsT1[5]=mxCreateDoubleScalar(nvec);
		prhsT1[6]=mxCreateDoubleScalar(T);
		// Generate the new vector using the matlab function "generator"
		mexCallMATLAB(nlhsT, plhsT1, nrhsT, prhsT1, "genParPowell");		
		// Recover the output
		//double *temp=mxGetPr(plhsT1[0]);
		//for(j=0;j<nvec;j++) mxSetCell(newV, j, mxGetCell(plhsT1[0],j));// newVec[j]=temp[j];
        newV=mxDuplicateArray(plhsT1[0]);
        //mexPrintf("#");
        //mxDestroyArray(plhsT1);
        //mxDestroyArray(prhsT1);
        
		
		//Generate and assign the input variable structures
        mxArray *prhsT2[6];
        mxArray *plhsT2[1];
		nlhsT=1, nrhsT=6;
		prhsT2[0]=test_stim; //stim
		prhsT2[1]=test_spikes; //spikes
		prhsT2[2]=newV;
		prhsT2[3]=nfilt; //nfilt
		prhsT2[4]=bins; //bins
		prhsT2[5]=mxCreateDoubleScalar(nvec);
		// Generate and recover the new negative information
		mexCallMATLAB(nlhsT, plhsT2, nrhsT, prhsT2, "findNegInfoMEX");
		testEnergy = mxGetScalar(plhsT2[0]);
        //mxDestroyArray(prhsT2);
        //mxDestroyArray(plhsT2);

		mxArray *prhsT8[6];
		mxArray *plhsT8[1];
		nlhsT=1, nrhsT=6;
		prhsT8[0]=stim; //stim
		prhsT8[1]=spikes; //spikes
		prhsT8[2]=newV;
		prhsT8[3]=nfilt; //nfilt
		prhsT8[4]=bins; //bins
		prhsT8[5]=mxCreateDoubleScalar(nvec);
		// Generate and recover the new negative information
		mexCallMATLAB(nlhsT, plhsT8, nrhsT, prhsT8, "findNegInfoMEX");
		newEnergy = mxGetScalar(plhsT8[0]);
        //mxDestroyArray(prhsT8);
        //mxDestroyArray(plhsT8);

		//mexPrintf("New Info: %5.5f\n",newEnergy);
		//mexPrintf("check1\n");
		// If the new energy is better than our best required...
		if(newEnergy < minEnergy)
		{
			mexPrintf("Better than maximum requirements, breaking\n\n");
			plhs[0] = mxDuplicateArray(bestV); //outVec[j]=newVec[j];
			//plhs[1] = mxCreateDoubleScalar(newEnergy);
			break;
		}
		
		// If our new energy is better than the best so far, update!
		else if(newEnergy < bestEnergy)
		{
            //mexPrintf("New Best\n");
			mexPrintf("%u: \t Updating... \n New Information: %7.8f,\t Current Temperature: %1.3e.\n",i, -newEnergy,T);
            FILE *out = fopen("anneal.txt", "a+");
			fprintf(out, "%f\n",newEnergy);
			fclose(out);
			//for(j=0;j<nvec;j++) mxSetCell(parentV, j, mxGetCell(newV,j));//parentVec[j]=newVec[j];
            parentV=mxDuplicateArray(newV);
			oldEnergy=newEnergy;
			success++;
			consec=0;
			//for(j=0;j<nvec;j++) mxSetCell(bestV, j, mxGetCell(newV,j));//bestVec[j]=newVec[j];
            bestV=mxDuplicateArray(newV);
			bestEnergy=newEnergy;
			plhs[0] = mxDuplicateArray(bestV);
			plhs[1]= mxCreateDoubleScalar(T);
		}
		
		// If our new energy is really bad, try a linear combination with the best vector
		else if(newEnergy-oldEnergy > 1000*T)
		{
			//mexPrintf("chooseVec\n");
			//mexPrintf("%u: \t Trying a scaled random vector offset...\n",i);
			nlhsT=1, nrhsT=8;

			//Assign the input variables for the function call
			mxArray *prhsT3[8];
			mxArray *plhsT3[1];
		    prhsT3[0]=newV;
		    prhsT3[1]=mxCreateDoubleScalar(1/T);
		    prhsT3[2]=stim; //stim
		    prhsT3[3]=spikes; //spikes
		    prhsT3[4]=nfilt; //nfilt
		    prhsT3[5]=bins; //bins
		    prhsT3[6]=mxCreateDoubleScalar(nvec);
	    		prhsT3[7]=mxCreateDoubleScalar(T);
		    //Call the matlab function "chooseVec.m"
		    mexCallMATLAB(nlhsT, plhsT3, nrhsT, prhsT3, "parPowellJump");
		    // Recover the output from the function call
		    //double *temp=mxGetPr(plhsT5[0]);
		    //for(j=0; j<nvec; j++) mxSetCell(parentV, j, mxGetCell(plhsT5[0],j));//parentVec[j]=temp[j];
            parentV=mxDuplicateArray(plhsT3[0]);
			//mxDestroyArray(plhsT3);
			//mxDestroyArray(prhsT3);

			//Generate and assign the input variable structures
            mxArray *prhsT4[6];
            mxArray *plhsT4[1];
			nlhsT=1, nrhsT=6;			
			prhsT4[0]=stim; //stim
			prhsT4[1]=spikes; //spikes
			prhsT4[2]=parentV;
			prhsT4[3]=nfilt; //nfilt
			prhsT4[4]=bins; //bins
			prhsT4[5]=mxCreateDoubleScalar(nvec);
			// Generate and recover the new negative information
			mexCallMATLAB(nlhsT, plhsT4, nrhsT, prhsT4, "findNegInfoMEX");
			oldEnergy = mxGetScalar(plhsT4[0]);
            //mxDestroyArray(prhsT4);
            //mxDestroyArray(plhsT4);

			consec++;
		}
		
		//  If the new energy is almost the same as the old energy, try shifting by a jump proportional to the current temperature
		else if(2*fabs(newEnergy-oldEnergy) <= (5.0e-6)*(fabs(newEnergy)+fabs(oldEnergy)))
		{
            //mexPrintf("jump\n");

                //mexPrintf("%u: delta-D_kl within tolerance, adding jitter.  Temperature: %1.3e\n",i,T);

			    // JUMP FUNCTION
			    nlhsT=1, nrhsT=8;
			    mxArray *prhsT5[8];
			    mxArray *plhsT5[1];
			    //Assign the input variables for the function call
			    prhsT5[0]=newV;
			    prhsT5[1]=mxCreateDoubleScalar(T);
			    prhsT5[2]=stim; //stim
			    prhsT5[3]=spikes; //spikes
			    prhsT5[4]=nfilt; //nfilt
			    prhsT5[5]=bins; //bins
			    prhsT5[6]=mxCreateDoubleScalar(nvec);
        	    		prhsT5[7]=mxCreateDoubleScalar(T);
			    //Call the matlab function "chooseVec.m"
			    mexCallMATLAB(nlhsT, plhsT5, nrhsT, prhsT5, "parPowellJump");
			    // Recover the output from the function call
			    //double *temp=mxGetPr(plhsT5[0]);
			    //for(j=0; j<nvec; j++) mxSetCell(parentV, j, mxGetCell(plhsT5[0],j));//parentVec[j]=temp[j];
                parentV=mxDuplicateArray(plhsT5[0]);
                //mxDestroyArray(prhsT5);
                //mxDestroyArray(plhsT5);
			
			    // ENERGY FUNCTION
			    //Generate and assign the input variable structures
			    nlhsT=1, nrhsT=6;
			    mxArray *prhsT6[6];
			    mxArray *plhsT6[1];
			    prhsT6[0]=stim; //stim
			    prhsT6[1]=spikes; //spikes
			    prhsT6[2]=parentV;
			    prhsT6[3]=nfilt; //nfilt
			    prhsT6[4]=bins; //bins
			    prhsT6[5]=mxCreateDoubleScalar(nvec);
			    // Generate and recover the new negative information
			    mexCallMATLAB(nlhsT, plhsT6, nrhsT, prhsT6, "findNegInfoMEX");
			    oldEnergy = mxGetScalar(plhsT6[0]);
                //mxDestroyArray(prhsT6);
                //mxDestroyArray(plhsT6);

			
			    // Heat up the temperature, if the heating coefficient is greater than one
			    if(heatCoef > 1.0)
			    {
				    T*=heatCoef;
				    if(T>maxTemp) T=maxTemp;
				    plhs[1] = mxCreateDoubleScalar(T);
			    }
			
			consec++;
		}
		// If nothing else held true...
		else
		{
            //mexPrintf("anneal\n");
			// Generate a random number via MATLAB
            mxArray *prhsT7[0], *plhsT7[1];
			nrhsT=0;
			nlhsT=1;
			mexCallMATLAB(nlhsT, plhsT7, nrhsT, prhsT7, "rand");
			double r=mxGetScalar(plhsT7[0]);
			//mxDestroyArray(prhsT7);
			//mxDestroyArray(plhsT7);
			//mexPrintf("Random NUmber: %4.4f\n",r);
			//If the random number is less than the exponential, keep trying along that gradient
			if(r < exp( (oldEnergy- newEnergy) / T))
			{
			    //for(j=0;j<nvec;j++) mxSetCell(parentV, j, mxGetCell(newV,j));//parentVec[j]=newVec[j];
                parentV=mxDuplicateArray(newV);
				oldEnergy=newEnergy;
				success++;
				consec=0;
			}
			
			// Otherwise, can it and try again
			else
			{
				consec++;
				//mexPrintf("%u: Settling... Temp: %1.4e\n",i,T);
			}



		}
	}
	// Assign the final output values
	//mexPrintf("\n Done! \n");

		
	//Free Memory

	mxDestroyArray(stim);
	mxDestroyArray(spikes);
	mxDestroyArray(nfilt);
	mxDestroyArray(bins);
	//mxDestroyArray(bestV);	
	mxDestroyArray(currentV);	
	mxDestroyArray(parentV);	
	mxDestroyArray(newV);
	mxDestroyArray(test_stim);
	mxDestroyArray(test_spikes);

}
