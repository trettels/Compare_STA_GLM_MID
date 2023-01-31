#include "mex.h"
#include <stdio.h>
#include <math.h>

/*
 *  linmin.c
 *  
 *
 *  Created by Sean Trettel on 6/30/10.
 *  Fairhall Lab, Department of Physiology and Biophysics, University of Washington
 *
 *	USAGE:
 *		[xmin]=linmin(vec, dvec, eps, maxIter, stim, spikes, nfilt, bins);
 *
 *	INPUT:
 *		vec = The filter-vector being examined
 *		dvec = the information gradient for vec
 *		eps = an epsilon value for tolerance creation
 *		maxIter = The maximum iterations to allow the Brent algorithm to run
 *		stim = the stimulus vector for passing on to findNegInfo
 *		spikes = the spike train vector for passing on to findNegInfo
 *		nfilt = the length of vec;
 *		bins = the number of bins to histogram over, passed on to findNegInfo
 *
 *	OUTPUT:
 *		xmin: the coefficient by which to multiply the gradient in [v'=v+xmin*dv] in order to minimize the negative information of v'.
 */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//mexPrintf("\n \t \t \t Entered Linear Minimization...\n");
	
	
	if(nrhs!=10) mexErrMsgTxt("Wrong number of input arguments"); //Error message thrown if there are an incorrect number of input arguments
	if(nlhs!=1) mexErrMsgTxt("Wrong number of output arguments"); //Error message thrown if there are an incorrect number of output arguments  
	
	
	// Initialize varaibles
	mxArray *v = mxDuplicateArray(prhs[0]);
    //mexPrintf(mxGetClassName(v));
    //mexPrintf("\n");
	mxArray *dv = mxDuplicateArray(prhs[1]);
    //mexPrintf(mxGetClassName(dv));
    //mexPrintf("\n");
	mxArray *testV = mxDuplicateArray(prhs[0]);
   // mexPrintf(mxGetClassName(testV));
    //mexPrintf("\n");
	//double *vec=mxGetPr(v);
	//double *dvec=mxGetPr(dv);
	//double *testVec = mxGetPr(testV);
	double eps = mxGetScalar(prhs[2]);
	double maxIter=mxGetScalar(prhs[3]);
	mxArray *stim = mxDuplicateArray(prhs[4]);
	mxArray *spikes = mxDuplicateArray(prhs[5]);
	mxArray *nfilt = mxDuplicateArray(prhs[6]);
	mxArray *bins = mxDuplicateArray(prhs[7]);
	int length = (int)(mxGetScalar(prhs[6]));
    int nvec=(int)(mxGetScalar(prhs[8]));
    int ihcurr= (int)(mxGetScalar(prhs[9]));

    

    //int HCURRENT = (int)(mxGetScalar(mexGetVariable("global","HCURRENT")));

	plhs[0]=mxCreateNumericMatrix(1,nvec,mxDOUBLE_CLASS,mxREAL);
    double *lambda = mxGetPr(plhs[0]);
	
	double x, a, b, fx, fa, fb, xmin, fmin, xminTemp, fminTemp;
	//mxArray *output1[1], *output2[1], *output3[1], *output4[1], *output5[1], *output6[1], *output7[1], *output8[1];
	//mxArray *input1[6], *input2[6], *input3[6], *input4[6], *input5[6], *input6[6], *input7[6], *input8[6];
	int i, j, q;
	double GOLD =  /* (3.0+sqrt(5.0))/2.0; */  1.618034;
	long maxDil = 300;
    //mexPrintf("Finished Initiallization, starting braketing\n");
    for(q=0;q<nvec;q++)
    {
        mxArray *vTemp=mxDuplicateArray(mxGetCell(v,q));
        //mexPrintf(mxGetClassName(vTemp));
        //mexPrintf("\n");
        mxArray *dvTemp = mxDuplicateArray(mxGetCell(dv,q));
        //mexPrintf(mxGetClassName(dvTemp));
        //mexPrintf("\n");
        mxArray *testVTemp=mxDuplicateArray(mxGetCell(testV,q));
       // mexPrintf(mxGetClassName(testVTemp));
        //mexPrintf("\n");
       // mexPrintf("reference vectors\n");
        double *vec = mxGetPr(vTemp);
        double *dvec = mxGetPr(dvTemp);
        double *testVec = mxGetPr(testVTemp);
        	// Bracketing Algorithm
        	a=0.0;
	    int onevlen=(int)(round(length/2));
        	double normV=0.0;
        	for(j=0; j< length; j++) normV += vec[j]*vec[j];
        	normV=sqrt(normV);
        	
        	double normDV=0.0;
        	for(j=0; j< length; j++) normDV += dvec[j]*dvec[j];
        	normDV=sqrt(normDV);
        	
        	x=-2.0;//*normV/(normDV+mxGetEps());
        	//mexPrintf("\n starting midpoint: %7.5f\n",x);
        //	mexPrintf("checkpoint1\n");
        	//Find the negative information for the left combination...
        	for(j=0; j<length; j++) testVec[j] = vec[j]+a*dvec[j];
            mxSetCell(testV, q, testVTemp);
        mxArray *output1[1], *input1[7];
        	int nlhsT=1, nrhsT=7;
        	input1[0]=stim; //stim
        	input1[1]=spikes; //spikes
        	input1[2]=testV;
        	input1[3]=nfilt; //nfilt
        	input1[4]=bins; //bins
        input1[5]=mxCreateDoubleScalar(nvec);
        input1[6]=mxCreateDoubleScalar(ihcurr);
        	// Generate and recover the negative information
        	mexCallMATLAB(nlhsT, output1, nrhsT, input1, "findNegInfoMEX");
        	fa = mxGetScalar(output1[0]);

        	//mxDestroyArray(input1);
        //mxDestroyArray(output1);
        //mexPrintf("a\n");
        
        	
        	//Find the negative information for the middle combination...
        	for(j=0; j<length; j++) testVec[j] = vec[j]+x*dvec[j];
        mxSetCell(testV, q, testVTemp);
        mxArray *output2[1], *input2[7];
        	nlhsT=1, nrhsT=7;
        	input2[0]=stim; //stim
        	input2[1]=spikes; //spikes
        	input2[2]=testV;
        	input2[3]=nfilt; //nfilt
        	input2[4]=bins; //bins
        input2[5]=mxCreateDoubleScalar(nvec);
        input2[6]=mxCreateDoubleScalar(ihcurr);
        	// Generate and recover the negative information
        	mexCallMATLAB(nlhsT, output2, nrhsT, input2, "findNegInfoMEX");
        	fx = mxGetScalar(output2[0]);

        //mxDestroyArray(input2);
        //mxDestroyArray(output2);
        //mexPrintf("fx\n");
        	
        	b=x+GOLD*(x-(a)); // =x;
        	//Find the negative information for the middle combination...
        	for(j=0; j<length; j++) testVec[j] = vec[j]+b*dvec[j];
        	mxSetCell(testV, q, testVTemp);
        mxArray *input8[7], *output8[1];
        	nlhsT=1, nrhsT=7;
        	input8[0]=stim; //stim
        	input8[1]=spikes; //spikes
        	input8[2]=testV;
        	input8[3]=nfilt; //nfilt
        	input8[4]=bins; //bins
        input8[5]=mxCreateDoubleScalar(nvec);
        input8[6]=mxCreateDoubleScalar(ihcurr);
        	// Generate and recover the negative information
        	mexCallMATLAB(nlhsT, output8, nrhsT, input8, "findNegInfoMEX");
        	fb = mxGetScalar(output8[0]); // =fx;
        	//mexPrintf("checkpoint2\n");

        //mxDestroyArray(input8);
        //mxDestroyArray(output8);
        	
        	if(fb < fa)
        	{
        		//mexPrintf("\n fb < fa\n");
        		b  *= 1+GOLD; 
        		//Find the negative information for the current combination...
        		for(j=0; j<length; j++) testVec[j] = vec[j]+b*dvec[j];
        		nlhsT=1, nrhsT=7;
        		mxArray *input3[7], *output3[1];
            mxSetCell(testV, q, testVTemp);
        		input3[0]=stim; //stim
        		input3[1]=spikes; //spikes
        		input3[2]=testV;
        		input3[3]=nfilt; //nfilt
        		input3[4]=bins; //bins
            input3[5]=mxCreateDoubleScalar(nvec);
            input3[6]=mxCreateDoubleScalar(ihcurr);
        		// Generate and recover the negative information
        		mexCallMATLAB(nlhsT, output3, nrhsT, input3, "findNegInfoMEX");
        		fb = mxGetScalar(output3[0]);

            //mxDestroyArray(input3);
            //mxDestroyArray(output3);
        
        		int iter=1;
        		
        		double min;
        		if(fa < fb) min=fa;
        		else min=fb;
        		
        		while(fx > min  && iter < maxDil)
        		{
        			x=b;
        			fx=fb;
        			b  *= 1+GOLD; 
        			//Find the negative information for the current combination...
        			for(j=0; j<length; j++) testVec[j] = vec[j]+b*dvec[j];
        			nlhsT=1, nrhsT=7;
        			mxArray *input4[7], *output4[1];
                mxSetCell(testV, q, testVTemp);
        			input4[0]=stim; //stim
        			input4[1]=spikes; //spikes
        			input4[2]=testV;
        			input4[3]=nfilt; //nfilt
        			input4[4]=bins; //bins
                input4[5]=mxCreateDoubleScalar(nvec);
                input4[6]=mxCreateDoubleScalar(ihcurr);
        			// Generate and recover the negative information
        			mexCallMATLAB(nlhsT, output4, nrhsT, input4, "findNegInfoMEX");
        			fb = mxGetScalar(output4[0]);

                //mxDestroyArray(input4);
                //mxDestroyArray(output4);
                
        			if(fa < fb) min=fa;
        			else min=fb;
        			
        			iter++;
        		}
        
        	}
        	else
        	{
        		//mexPrintf("\n fb > fa\n");
        		x /= 1+GOLD; 
        		//Find the negative information for the current combination...
        		for(j=0; j<length; j++) testVec[j] = vec[j]+x*dvec[j];
        		nlhsT=1, nrhsT=7;
        		mxArray *input5[7], *output5[1];
            mxSetCell(testV, q, testVTemp);
        		input5[0]=stim; //stim
        		input5[1]=spikes; //spikes
        		input5[2]=testV;
        		input5[3]=nfilt; //nfilt
        		input5[4]=bins; //bins
            input5[5]=mxCreateDoubleScalar(nvec);
            input5[6]=mxCreateDoubleScalar(ihcurr);
        		// Generate and recover the negative information
        		mexCallMATLAB(nlhsT, output5, nrhsT, input5, "findNegInfoMEX");
        		fx = mxGetScalar(output5[0]);

            //mxDestroyArray(input5);
            //mxDestroyArray(output5);
        
        		int iter=1;
        		double min;
        		if(fa <= fb) min=fa;
        		else min=fb;
        		while(fx>min && iter < maxDil)
        		{
        			b=x;
        			fb=fx;
        			x  /= 1+GOLD;
        			//Find the negative information for the current combination...
        			for(j=0; j<length; j++) testVec[j] = vec[j]+x*dvec[j];
        			nlhsT=1, nrhsT=7;
        			mxArray *input6[7], *output6[1];
                mxSetCell(testV, q, testVTemp);
                input6[0]=stim; //stim
        			input6[1]=spikes; //spikes
        			input6[2]=testV;
        			input6[3]=nfilt; //nfilt
        			input6[4]=bins; //bins
                input6[5]=mxCreateDoubleScalar(nvec);
                input6[6]=mxCreateDoubleScalar(ihcurr);
        			// Generate and recover the negative information
        			mexCallMATLAB(nlhsT, output6, nrhsT, input6, "findNegInfoMEX");
        			fx = mxGetScalar(output6[0]);

                //mxDestroyArray(input6);
                //mxDestroyArray(output6);
        
        			if(fa <= fb) min=fa;
        			else min=fb;
        			
        			iter++;
        		}
        
        
        	}
        	//mexPrintf("\n Finished Bracketing, entering Brent's Algorithm\n");
        	//mexPrintf("Finished Bracket, Starting Brent\n");
        	//Brent's Algorithm
        	double e=0.0;
        	fmin = fx;
        	fminTemp = fx;
        	xmin = x;
        	xminTemp = x;
        	double Zeps = 1.0e-10;
        
        	
        	for(i=0; i<maxIter; i++)
        	{
        		double mid=0.5*(a+b);
        		double tol1 = eps*fabs(xmin)+Zeps;
        		double tol2=2.0*tol1;
        		double u, d, fu;
        		
        		if(fabs(xmin-mid) <= tol2-0.5*(b-a)) // We're done!
        		{
        			lambda[q]=xmin;
        			//mexPrintf("\n Exit Linear Minimization early\n ");
        			break;
        		}
        		else if(fabs(e) <= tol1)
        		{
        			if(xmin >= mid) e=a-xmin;
        			else e=b-xmin;
        			d=e/GOLD;
        		}
        		else //Parabolic test-fit
        		{
        			double r = (xmin-xminTemp)*(fmin-fx);
        			double q = (xmin-x)*(fmin-fminTemp);
        			double p=(xmin-x)*q-(xmin-xminTemp)*r;
        			q=2.0*(q-r);
        			
        			if(q>0.0) p=-p;
        			q=fabs(q);
        			double eTemp=e;
        			e=d;
        			
        			if(fabs(p) >= fabs(0.5*q*eTemp) || p<=q*(a-xmin) || p>=q*(b-xmin))
        			{
        				if(xmin>x) e=a-xmin;
        				else e=b-xmin;
        				d=e/GOLD;
        			}
        			else
        			{
        				d=p/q;
        				u=xmin+d;
        				if(u-a < tol2 || b-u <tol2 && (x-xmin) != 0.0) d=tol1*(x-xmin)/fabs(x-xmin);
        				else d=0.0;
        			}
        		}
        		if(fabs(d) >= tol1) u=xmin+d;
        		else if(d!=0) u=xmin+tol1*d/fabs(d);
        		else u=xmin;
        		
        		//Find the negative information for the current combination...
        		for(j=0; j<length; j++) testVec[j] = vec[j]+u*dvec[j];
        		nlhsT=1, nrhsT=7;
        		mxArray *input7[7], *output7[1];
            mxSetCell(testV, q, testVTemp);
        		input7[0]=stim; //stim
        		input7[1]=spikes; //spikes
        		input7[2]=testV;
        		input7[3]=nfilt; //nfilt
        		input7[4]=bins; //bins
            input7[5]=mxCreateDoubleScalar(nvec);
            input7[6]=mxCreateDoubleScalar(ihcurr);
        		// Generate and recover the negative information
        		mexCallMATLAB(nlhsT, output7, nrhsT, input7, "findNegInfoMEX");
        		fu = mxGetScalar(output7[0]);

            //mxDestroyArray(input7);
            //mxDestroyArray(output7);
        		
        		if(fu <= fmin)
        		{
        			if(u>=xmin) a=xmin;
        			else b=xmin;
        			
        			x=xminTemp;
        			fx=fminTemp;
        			xminTemp =xmin;
        			fminTemp = fmin;
        			xmin=u;
        			fmin=fu;
        		}
        		else
        		{
        			if(u < xmin) a=u;
        			else b=u;
        			
        			if(fu <=fminTemp)
        			{
        				x=fminTemp;
        				fx=fminTemp;
        				xminTemp=u;
        				fminTemp=fu;
        			}
        			else if(fu<=fx || x==xmin || xminTemp == x)
        			{
        				x=u;
        				fx=fu;
        			}
        		}
        	}
    
    	lambda[q] = xmin;
    	
    }
    //mexPrintf("exit linmin");
    
        	//Free Memory
    	mxDestroyArray(v);
    	mxDestroyArray(dv);
    	mxDestroyArray(testV);
    	mxDestroyArray(stim);
    	mxDestroyArray(spikes);
    	mxDestroyArray(nfilt);
    	mxDestroyArray(bins);
     
}  
