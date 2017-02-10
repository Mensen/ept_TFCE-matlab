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

 *
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

void tfce_thread(double *inData, double *outData, double *ChN, double *thresh, const int *dims, const int *dims2)
{
   double valToAdd;
   int i, t, f, ti, tt, tf, maxt, mint, maxf, minf, temp, growingInd, growingCur, ChCurr, idx, p=0, n=0;
   int numVoxels = dims[0] * dims[1] * dims[2];
   char* flagUsed;
   /*#-short* growing; */
   short* grow_i;
   short* grow_t;
   short* grow_f;
   
   
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
				/* temp is the current point in the grid */
				temp = (f*dims[1]*dims[0]) + (t*dims[0]) + i;
				
				/* Check if this point has been seen (flagUsed) and whether its over the threshold */
				if (!flagUsed[temp] && inData[temp] >= thresh[0])
				{
				    /* mexPrintf("\nPositive Values... \n");
					/* mexPrintf("Looking at %f \n", inData[temp]);
					/* mexPrintf("t is %d \n", t);
					/* mexPrintf("f is %d \n", f);
					/* mexEvalString("drawnow;"); */
								
					/* make the flagUsed so that algorithm doesn't visit this point again */
					flagUsed[temp] = 1;
					growingInd = 1;
					growingCur = 0;
					
					/* Define the current coordinates to continue with */
					grow_i[0]=i;
					grow_t[0]=t;
					grow_f[0]=f;
					p++;
					
					while (growingCur < growingInd)
					{
					   /*This just limits not to overrun borders <-- in our case the zero padding
					   /*And creates a 3x3 windows of scanning for valid points (one point around the current point)
					   /*E.g. In 2 dimensions... if the current point is 2,j then maxi = Min of dimension size or 2+2
					   /*thus maxi = 4... mini = max of 0 or 2-1... so 1...
					   /*so in the next set of loops ti will look at i=1,i=2,and i=3 (not 4 because its set to < maxi */
					   maxt = MIN(dims[1], grow_t[growingCur] + 2);
					   maxf = MIN(dims[2], grow_f[growingCur] + 2);
					   mint = MAX(0, grow_t[growingCur] - 1);
					   minf = MAX(0, grow_f[growingCur] - 1);
						
						/* start of the smaller scanning window */
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
									
									if (!flagUsed[temp] && inData[temp] >= thresh[0])
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
					   /* GrowingCur increases one and thus looks at the next point that was found in the previous loop */
					   growingCur++;

					}
					/* Reset back to 0 so that when the next "while loop" runs it adds the value to each point found */
					growingCur = 0;
								
					valToAdd = p;
					
					while (growingCur < growingInd)
					{
						outData[(grow_f[growingCur]*dims[1]*dims[0]) + (grow_t[growingCur]*dims[0]) + grow_i[growingCur]] += valToAdd;

					   growingCur++;
					}
				}
				
				temp = (f*dims[1]*dims[0]) + (t*dims[0]) + i;
				if (!flagUsed[temp] && -inData[temp] >= thresh[0])
				{
					flagUsed[temp] = 1;
					growingInd = 1;
					growingCur = 0;
					
					grow_i[0]=i;
					grow_t[0]=t;
					grow_f[0]=f;
					n++;

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
									if (!flagUsed[temp] && -inData[temp] >= thresh[0])
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
								
					valToAdd = n;
					
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

/* Declarations */
double *inData, *outData, *ChN, *thresh; 
int ndim, ndim2;
const int *dims, *dims2;

/* check inputs */
if (nrhs!=3)
  mexErrMsgTxt("3 inputs required. (Data, ChN, Thresh");
else if (nlhs>1)
  mexErrMsgTxt("Too many output arguments.");
  
if (!mxIsDouble(prhs[0]))
	mexErrMsgTxt("First argument must be double.");

/* get input inDatage */
inData = (double*)mxGetPr(prhs[0]);
ChN    = (double*)mxGetPr(prhs[1]);
thresh = (double*)mxGetPr(prhs[2]);

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

tfce_thread(inData, outData, ChN, thresh, dims, dims2); 

return;
}

