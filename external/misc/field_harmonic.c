/* 
Copyright © 2011 by Université catholique de Louvain.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files, to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

REGGUI SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/*#include <omp.h>*/
#include <math.h>
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Inputs */
	float	*datf;	/* field */
	float	*spaf;	/* spacings */

	/* Outputs */
	mxArray *harmonic;
	float	*dath;	/* harmonic energy */

	/* Size vector */
	const int	*sizf;

	/* Variables */
	long	i, j, k;		/* Indices */
	long	ubif, ubjf, ubkf;	/* Upper bounds */
	long	osif, osjf, oskf;	/* Offsets */
	long	osxf, osyf, oszf;	/* Offsets */
	long	indf;			/* Index */
	float	diff;
	float	mat9[9];		/* 3-by-3 Jacobian matrix */
	float	*ptxf, *ptyf, *ptzf;	/* pointers to current entry in the three components of the field */

	/* Inputs */

	/* Check the number of input and output arguments. */
	if (nrhs != 2) mexErrMsgTxt("This function takes two arguments.");

	/* We won't deal with empty arrays. */
	for (i = 0; i < nrhs; i++)
	{	if (mxIsEmpty(prhs[i])) mexErrMsgTxt("Arguments may not be empty.");
	}

	/* Check the forspaf of the input arguments. */
	for (i = 0; i < nrhs; i++)
	{	if (!mxIsNumeric(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsSingle(prhs[i]) || mxIsComplex(prhs[i]))
		{	mexErrMsgTxt("All input arguments must be full numeric arrays of real numbers, stored as singles.");
		}
	}

	/* Get the field size. Trailing singleton dimensions are ignored. */
	sizf = mxGetDimensions(prhs[0]);
	for (i = mxGetNumberOfDimensions(prhs[0]); i > 0; i--) if (sizf[i-1] > 1) break;

	/* Verify that there are exactly 4 dimensions. */
	if (4 != i) mexErrMsgTxt("The first argument (deformation field) must have 4 dimensions.");

	/* Get the spacing size */
	if (3 != mxGetM(prhs[1])*mxGetN(prhs[1])) mexErrMsgTxt("The second argument (spacings) must have 3 elements");
	spaf = (float*)mxGetPr(prhs[1]);

	/* Output */

	/* Create the output array. */
	harmonic = mxCreateNumericArray(3, sizf, mxSINGLE_CLASS, mxREAL);

	/* Get pointer to the field. */
	datf = (float*)mxGetPr(prhs[0]);

	/* Get pointer to the harmonic. */
	dath = (float*)mxGetPr(harmonic);

	/* OpenMP */
#ifdef _OPENMP
	omp_set_num_threads(4);
	/* omp_set_dynamic(1); */
#endif

	ubif = (long)sizf[0] - 1;
	ubjf = (long)sizf[1] - 1;
	ubkf = (long)sizf[2] - 1;

	osif = 1;
	osjf = (long)sizf[0] * osif;
	oskf = (long)sizf[1] * osjf;

	osxf = 0;
	osyf = (long)sizf[2] * oskf;
	oszf = 2 * osyf;

	#pragma omp parallel for private(i,j,indf,diff,mat9)
	for (k = 0; k <= ubkf; k++)
	{
		for (j = 0; j <= ubjf; j++)
		{
			for (i = 0; i <= ubif; i++)
			{	/* index */
				indf = i*osif + j*osjf + k*oskf;

				/* pointers to the three components */
				ptxf = datf + indf + osxf;
				ptyf = datf + indf + osyf;
				ptzf = datf + indf + oszf;

				/* compute Jacobian */

				/* partial derivative w.r.t. x */
				diff = 0.0;

				if (ubif>i)
				{	mat9[0] = *(ptxf+osif);
					mat9[1] = *(ptyf+osif);
					mat9[2] = *(ptzf+osif);
					diff += spaf[0];
				}
				else
				{	mat9[0] = *ptxf;
					mat9[1] = *ptyf;
					mat9[2] = *ptzf;
				}

				if (0<i)
				{	mat9[0] -= *(ptxf-osif);
					mat9[1] -= *(ptyf-osif);
					mat9[2] -= *(ptzf-osif);
					diff += spaf[0];
				}
				else
				{	mat9[0] -= *ptxf;
					mat9[1] -= *ptyf;
					mat9[2] -= *ptzf;
				}

				if (diff)
				{	mat9[0] += diff;	/* diagonal term */
					mat9[0] /= diff;
					mat9[1] /= diff;
					mat9[2] /= diff;
				}

				/* partial derivative w.r.t. y */
				diff = 0.0;

				if (ubjf>j)
				{	mat9[3] = *(ptxf+osjf);
					mat9[4] = *(ptyf+osjf);
					mat9[5] = *(ptzf+osjf);
					diff += spaf[1];
				}
				else
				{	mat9[3] = *ptxf;
					mat9[4] = *ptyf;
					mat9[5] = *ptzf;
				}

				if (0<j)
				{	mat9[3] -= *(ptxf-osjf);
					mat9[4] -= *(ptyf-osjf);
					mat9[5] -= *(ptzf-osjf);
					diff += spaf[1];
				}
				else
				{	mat9[3] -= *ptxf;
					mat9[4] -= *ptyf;
					mat9[5] -= *ptzf;
				}

				if (diff)
				{	mat9[4] += diff;	/* diagonal term */
					mat9[3] /= diff;
					mat9[4] /= diff;
					mat9[5] /= diff;
				}

				/* partial derivative w.r.t. y */
				diff = 0.0;

				if (ubkf>k)
				{	mat9[6] = *(ptxf+oskf);
					mat9[7] = *(ptyf+oskf);
					mat9[8] = *(ptzf+oskf);
					diff += spaf[2];
				}
				else
				{	mat9[6] = *ptxf;
					mat9[7] = *ptyf;
					mat9[8] = *ptzf;
				}

				if (0<k)
				{	mat9[6] -= *(ptxf-oskf);
					mat9[7] -= *(ptyf-oskf);
					mat9[8] -= *(ptzf-oskf);
					diff += spaf[2];
				}
				else
				{	mat9[6] -= *ptxf;
					mat9[7] -= *ptyf;
					mat9[8] -= *ptzf;
				}

				if (diff)
				{	mat9[8] += diff;	/* diagonal term */
					mat9[6] /= diff;
					mat9[7] /= diff;
					mat9[8] /= diff;
				}

				/* determinant */
				dath[indf] = ( ((mat9[0]-1.0)*(mat9[0]-1.0)) + (mat9[1]*mat9[1]) + (mat9[2]*mat9[2]) + (mat9[3]*mat9[3]) + ((mat9[4]-1.0)*(mat9[4]-1.0)) + (mat9[5]*mat9[5]) + (mat9[6]*mat9[6]) + (mat9[7]*mat9[7]) + ((mat9[8]-1.0)*(mat9[8]-1.0)) ) /2.0;
			}
		}
	}

	/* Output the computed result. */
	plhs[0] = harmonic;
}


