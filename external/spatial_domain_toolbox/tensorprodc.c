#include "mex.h"

/*
 * TENSORPRODC
 * 
 * See the help text in tensorprodc.m for a description.
 * 
 * Author: Gunnar Farnebäck
 *         Medical Informatics
 *         Linköping University, Sweden
 *         gunnar@imt.liu.se
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int number_of_arrays = (nrhs - 2) / 2;
    double *index_sizes = mxGetPr(prhs[0]);
    int number_of_indices = mxGetM(prhs[0]) * mxGetN(prhs[0]);
    int output_ndims = mxGetM(prhs[1]) * mxGetN(prhs[1]);
    int *output_dims;
    double *output;
    int k;
    int *indices;
    int number_index_values = 1;
    int x;
    int max_array_dims = output_ndims;
    int *array_indices;
    int output_index;
    double p;
    int m;

    /* Create the result array.
     *
     * In case all indices are repeated the output is a scalar which
     * should logically be considered to be zero-dimensional, but that
     * doesn't go so well with Matlab's model of arrays. Instead we
     * have to special case it.
     */
    if (output_ndims == 0)
	plhs[0] = mxCreateDoubleScalar(0.0);
    else {
	output_dims = mxCalloc(output_ndims, sizeof(output_dims[0]));
	for (k = 0; k < output_ndims; k++) {
	    int n = (int) mxGetPr(prhs[1])[k] - 1;
	    output_dims[k] = (int) mxGetPr(prhs[0])[n];
	}
	plhs[0] = mxCreateNumericArray(output_ndims, output_dims,
				       mxDOUBLE_CLASS, mxREAL);
	mxFree(output_dims);
    }
    output = mxGetPr(plhs[0]);

    /* Allocate space to hold the index values and find out the total
     * number of index combinations to loop over.
     */
    indices = mxCalloc(number_of_indices, sizeof(indices[0]));
    for (k = 0; k < number_of_indices; k++)
	number_index_values *= index_sizes[k];

    /* Find the maximum number of indices for a single array. */
    for (m = 0; m < number_of_arrays; m++)
	if (mxGetNumberOfDimensions(prhs[m + 2]) > max_array_dims)
	    max_array_dims = mxGetNumberOfDimensions(prhs[m + 2]);

    /* Allocate space for the index values used by a single array.
     * This is shared between all arrays.
     */
    array_indices = mxCalloc(max_array_dims, sizeof(array_indices[0]));

    /* Loop over all combinations of index values. */
    for (x = 0; x < number_index_values; x++) {
	/* Compute index values for this iteration of the loop. */
	int y = x;
	for (k = 0; k < number_of_indices; k++) {
	    indices[k] = y % (int) index_sizes[k];
	    y /= index_sizes[k];
	}

	/* Find the linear subscript into the output array. Once more
	 * the scalar output case must be special-cased.
	 */
	if (output_ndims == 0)
	    output_index = 0;
	else {
	    for (k = 0; k < output_ndims; k++)
		array_indices[k] = indices[(int) mxGetPr(prhs[1])[k] - 1];
	    output_index = mxCalcSingleSubscript(plhs[0], output_ndims,
						 array_indices);
	}

	/* Compute the product of the array values for the current set
	 * of indices.
	 */
	p = 1.0;
	for (m = 0; m < number_of_arrays; m++) {
	    int n = mxGetM(prhs[2 + number_of_arrays + m]);
	    double *i = mxGetPr(prhs[2 + number_of_arrays + m]);
	    int index;
	    for (k = 0; k < n; k++)
		array_indices[k] = indices[(int) i[k] - 1];
	    index = mxCalcSingleSubscript(prhs[2 + m], n, array_indices);
	    p *= mxGetPr(prhs[2 + m])[index];
	}

	/* Sum this product to the result. */
	output[output_index] += p;
    }

    mxFree(indices);
    mxFree(array_indices);
}
