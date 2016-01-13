       /* This file is part of the program ept_TFCE.
	* ept_TFCE is free software: you can redistribute it and/or modify
	* it under the terms of the GNU General Public License as published by
	* the Free Software Foundation, either version 3 of the License, or
	* (at your option) any later version.
	* ept_TFCE is distributed in the hope that it will be useful,
	* but WITHOUT ANY WARRANTY; without even the implied warranty of
	* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	* GNU General Public License for more details.
	* You should have received a copy of the GNU General Public License
	* along with ept_TFCE.  If not, see <http://www.gnu.org/licenses/>. 
	*/

/*
 * Christian Gaser
 *
 * Pau Coma 09.12.2010
 *
 * -Independent grow_i, grow_j, grow_k variables to increase understandability
 *
 * Armand Mensen 
 * 20.12.2010
 * - Added the input which indicates channels neighbours
 * - Searches the 2D surface coordinates for neighbours
 *
 * 10.09.2013
 * - Searches over 3rd dimension now (Channel*Time*Frequency)
 *
 */

#include "math.h"
#include "mex.h"
#include <stdlib.h>

#ifdef OPENMP
#include "omp.h"
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

void tfce_thread(double *inData, double *outData, double *ChN, double thresh, double delta, const int *dims, const int *dims2, double *EH)
{
   double valToAdd;
   double E, H;
   int i, t, f, ti, tt, tf, maxt, mint, maxf, minf, temp, growingInd, growingCur, ChCurr, idx;
   int numVoxels = dims[0] * dims[1] * dims[2];
   char* flagUsed;
   short* grow_i;
   short* grow_t;
   short* grow_f;
   
   E = EH[0];
   H = EH[1];
   
   flagUsed = (char*)malloc(numVoxels*sizeof(char));
   grow_i  = (short*)malloc(numVoxels*sizeof(short));
   grow_t  = (short*)malloc(numVoxels*sizeof(short));
   grow_f  = (short*)malloc(numVoxels*sizeof(short));
   

    for (temp = 0; temp < numVoxels; ++temp) flagUsed[temp] = 0;
    
	for (f = 0; f < dims[2]; ++f)
	{
		for (t = 0; t < dims[1]; ++t)
		{
			for (i = 0; i < dims[0]; ++i)
			{
				temp = (f*dims[1]*dims[0]) + (t*dims[0]) + i;
				if (!flagUsed[temp] && inData[temp] >= thresh)
				{
					flagUsed[temp] = 1;
					growingInd = 1;
					growingCur = 0;
					
					grow_i[0]=i;
					grow_t[0]=t;
					grow_f[0]=f;

					while (growingCur < growingInd)
					{
					   maxt = MIN(dims[1], grow_t[growingCur] + 2);
					   maxf = MIN(dims[2], grow_f[growingCur] + 2);
					   mint = MAX(0, grow_t[growingCur] - 1);
					   minf = MAX(0, grow_f[growingCur] - 1);
						
						for (tf = minf; tf < maxf; ++tf)
						{
							for (tt = mint; tt < maxt; ++tt)
							{
								for (ti = 0; ti < dims2[1]; ++ti)
								{
									idx = (ti*dims2[0]) + grow_i[growingCur];
									ChCurr = ChN[idx];
									
									if (ChCurr == 0)
									{
										break;
									}
									
									ChCurr = ChCurr - 1;
									
									temp = tf*(dims[0]*dims[1]) + (tt*dims[0]) + ChCurr;
									
									if (!flagUsed[temp] && inData[temp] >= thresh)
									{
										flagUsed[temp] = 1;
										grow_i[growingInd] = ChCurr;
										grow_t[growingInd] = tt;
										grow_f[growingInd] = tf;
										growingInd++;
									}
								}	
							}
						}
					   growingCur++;

					}
					growingCur = 0;
								
					valToAdd = pow(growingInd, E) * pow(thresh, H) * delta;
					
					while (growingCur < growingInd)
					{
					   outData[(grow_f[growingCur]*dims[1]*dims[0]) + (grow_t[growingCur]*dims[0]) + grow_i[growingCur]] += valToAdd;
					   growingCur++;
					}
				}
				
/* Search for negative values... */

				if (!flagUsed[temp] && -inData[temp] >= thresh)
				{
					flagUsed[temp] = 1;
					growingInd = 1;
					growingCur = 0;
					
					grow_i[0]=i;
					grow_t[0]=t;
					grow_f[0]=f;

					while (growingCur < growingInd)
					{
					   maxt = MIN(dims[1], grow_t[growingCur] + 2);
					   maxf = MIN(dims[2], grow_f[growingCur] + 2);
					   mint = MAX(0, grow_t[growingCur] - 1);
					   minf = MAX(0, grow_f[growingCur] - 1);
						
						for (tf = minf; tf < maxf; ++tf)
						{
							for (tt = mint; tt < maxt; ++tt)
							{
								for (ti = 0; ti < dims2[1]; ++ti)
								{
									idx = (ti*dims2[0]) + grow_i[growingCur];
									ChCurr = ChN[idx];
									
									if (ChCurr == 0)
									{
										break;
									}
									
									ChCurr = ChCurr - 1;
									
									temp = tf*(dims[0]*dims[1]) + (tt*dims[0]) + ChCurr;
									
									if (!flagUsed[temp] && -inData[temp] >= thresh)
									{
										flagUsed[temp] = 1;
										grow_i[growingInd] = ChCurr;
										grow_t[growingInd] = tt;
										grow_f[growingInd] = tf;
										growingInd++;
									}
								}	
							}
						}
					   growingCur++;

					}
					growingCur = 0;
								
					valToAdd = pow(growingInd, E) * pow(thresh, H) * delta;
					
					while (growingCur < growingInd)
					{
					   outData[(grow_f[growingCur]*dims[1]*dims[0]) + (grow_t[growingCur]*dims[0]) + grow_i[growingCur]] -= valToAdd;
					   growingCur++;
					}
				}
			}
		}
	}
   
   free(flagUsed);
   
   free(grow_i);
   free(grow_t);
   free(grow_f);
}

void tfce(double *inData, double *outData, double *ChN, const int *dims, const int *dims2, double *EH)
{
   double fmax = 0.0, thresh0, thresh, delta;
   int i, numSteps = 50;
   int numVoxels = dims[0] * dims[1] * dims[2];
   
   for (i = 0; i < numVoxels; ++i)
   {
      if (inData[i] > fmax) fmax = inData[i];
	  else if (-inData[i] > fmax) fmax = -inData[i];
      outData[i] = 0.0;
   }
   
   delta = fmax/numSteps;
   thresh0 = delta/2.0;
   
  {
	#ifdef OPENMP
		/* this is faster than the default setting to the # of processors */
		omp_set_num_threads(numSteps);
		# pragma omp parallel for private(thresh) shared(outData) 
	#endif
    for (i = 0; i < numSteps; i++)
    {
       thresh = thresh0 + (double)i*delta;   
       tfce_thread(inData, outData, ChN, thresh, delta, dims, dims2, EH);
    }
  }   
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

/* Declarations */
double *inData, *outData, *ChN; 
double *EH;
int ndim, ndim2;
const int *dims, *dims2;

/* check inputs */
if (nrhs!=3)
  mexErrMsgTxt("3 inputs required. (Data, ChN, [E, H])");
else if (nlhs>1)
  mexErrMsgTxt("Too many output arguments.");
  
if (!mxIsDouble(prhs[0]))
	mexErrMsgTxt("First argument must be double.");

/* get input inDatage */
inData = (double*)mxGetPr(prhs[0]);
ChN    = (double*)mxGetPr(prhs[1]);
EH 	   = (double*)mxGetPr(prhs[2]);

ndim = mxGetNumberOfDimensions(prhs[0]);
ndim2 = mxGetNumberOfDimensions(prhs[1]);
if (ndim!=3){
	mexErrMsgTxt("Data input should be 3D: Channel*Frequency*Time");
}
else if (ndim2!=2){
	mexErrMsgTxt("Channel*Neighbour file should be 2D");
}
	
dims = mxGetDimensions(prhs[0]);
dims2 = mxGetDimensions(prhs[1]);

/*Allocate memory and assign output pointer*/
plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL);

/*Get a pointer to the data space in our newly allocated memory*/
outData = mxGetPr(plhs[0]);

tfce(inData, outData, ChN, dims, dims2, EH); 

return;
}

