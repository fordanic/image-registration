#include "mex.h"

#define DEBUG 0

/*
 * normconv.c
 *
 * Perform normalized convolution. Only real images and filters are
 * supported, but at arbitrary dimensionality. Perform normalized
 * convolution. Complex signals and basis functions are supported, at
 * arbitrary dimensionality. The method is described in chapter 3 of
 * Gunnar Farnebäck's thesis "Polynomial Expansion for Orientation and
 * Motion Estimation".
 *
 * Call:
 * result = normconv(signal, certainty, basis, applicability)
 * or
 * result = normconv(signal, certainty, basis, applicability,
 *                   region_of_interest)
 * 
 * signal             - signal values
 * certainty          - signal certainty (pointwise)
 * basis              - subspace basis
 * applicability      - applicability function for the basis
 * region_of_interest - where to compute the normalized convolution
 * result             - local basis coefficients
 * 
 * Formats: (Also see note 2 for special cases.)
 * signal is an n-dimensional array.
 * certainty is the same size as signal.
 * basis is an (n+1)-dimensional array. The first n dimensions correspond
 *         to the signal dimensions. The size of the last dimension gives
 *         the number of basis functions.
 * applicability is an n-dimensional array, with the same size as the first
 *         n dimensions of basis.
 * region_of_interest is either an nx2 matrix, where each row gives start
 *         and end indices for each dimension, or a sparse matrix of the
 *         same size as signal.
 * result is an (n+1)-dimensional array. The size of the first n dimensions
 *         are the same as for signal, unless region_of_interest is
 *         specified. The size of the last dimension equals the number of
 *         basis functions.
 * 
 * Note 1: Only double, nonsparse arrays are currently supported.
 * 
 * Note 2: Trailing singleton dimensions cannot exist in matlab 5.
 *         If there is a mismatch in the number of dimensions of the
 *         parameters it will be assumed that there are additional
 *         singleton dimensions, if that turns out to make sense.
 *         Particularly, in the case of a single basis function, basis
 *         will normally be n-dimensional instead of (n+1)-dimensional.
 *         The same goes for result.
 * 
 *         In the case of column vectors, the second, singleton, dimension
 *         will be ignored, although it is included by matlab 5.
 * 
 * Note 3: The special cases of 1-3 dimensions and/or a single basis
 *         function have been optimized for speed. To use the general
 *         algorithm also in these cases, add as a final parameter in
 *         the call the string 'general'. This will certainly be slower
 *         but can be used to verify that the optimized code works as
 *         intended. The results should not differ by more than
 *         numerical deviations. The general algorithm does not work with
 *         the sparse form of region_of_interest.
 * 
 * Note 4: The above mentioned sparse region of interest has not been
 *         implemented, as it turned out that matlab only supports
 *         sparse matrices in two dimensions.
 * 
 * Author: Gunnar Farnebäck
 *         Computer Vision Laboratory
 *         Linköping University, Sweden
 *         gf@isy.liu.se
 */


#define NO_ROI   0
#define MASK_ROI 1
#define BOX_ROI  2

/* Set up matrices needed for computation of output certainty. */
static void
init_cout_matrices(const double *basis_r, const double *basis_i,
		   const double *applicability,
		   int basis_size, int modeldimprod,
		   int *cout_num_arguments, mxArray **cout_arguments,
		   mxArray *G_matrix, mxArray *G0_matrix,
		   mxArray *h_vector, mxArray *cout_data_array,
		   int is_real)
{
  int k1, k2;
  double p_r, p_i;
  double ip_r, ip_i;
  int model_index;
  double *G0_r;
  double *G0_i;
  
  if (is_real)
  {
    G0_matrix = mxCreateDoubleMatrix(basis_size, basis_size, mxREAL);
    G0_r = mxGetPr(G0_matrix);
    G0_i = NULL;
  }
  else
  {
    G0_matrix = mxCreateDoubleMatrix(basis_size, basis_size, mxCOMPLEX);
    G0_r = mxGetPr(G0_matrix);
    G0_i = mxGetPi(G0_matrix);
  }
  
  /* Double loop over the basis functions to compute inner
   * products between basis functions for certainty identically one.
   */
  for (k1 = 0; k1 < basis_size; k1++)
  {
    for (k2 = k1; k2 < basis_size; k2++)
    {
      /* Reset the accumulated inner product. */
      ip_r = 0.0;
      ip_i = 0.0;
      
      /* Loop over the dimensions for the basis functions. */
      for (model_index = 0;
	   model_index < modeldimprod;
	   model_index++)
      {
	if (is_real)
	{
	  p_r = basis_r[model_index + modeldimprod * k2];
	  p_r *= applicability[model_index];
	  p_r *= basis_r[model_index + modeldimprod * k1];
	  ip_r += p_r;
	}
	else
	{
	  double b_r, b_i;
	  p_r = basis_r[model_index + modeldimprod * k2];
	  p_i = basis_i[model_index + modeldimprod * k2];
	  p_r *= applicability[model_index];
	  p_i *= applicability[model_index];
	  b_r = basis_r[model_index + modeldimprod * k1];
	  b_i = basis_i[model_index + modeldimprod * k1];
	  ip_r += p_r * b_r + p_i * b_i;
	  ip_i += p_i * b_r - p_r * b_i;
	}
      }
      
      if (is_real)
      {
	G0_r[k1 + k2 * basis_size] = ip_r;
	G0_r[k2 + k1 * basis_size] = ip_r;
      }
      else
      {
	G0_r[k1 + k2 * basis_size] = ip_r;
	G0_i[k1 + k2 * basis_size] = ip_i;
	G0_r[k2 + k1 * basis_size] = ip_r;
	G0_i[k2 + k1 * basis_size] = -ip_i;
      }
    }
  }
  
  *cout_num_arguments = 4;
  cout_arguments[0] = G_matrix;
  cout_arguments[1] = G0_matrix;
  cout_arguments[2] = h_vector;
  if (cout_data_array)
  {
    cout_arguments[4] = cout_data_array;
    *cout_num_arguments = 5;
  }
}


/* One signal dimension, multiple basis functions. */
static void
normconv1(const double *signal_r,
	  const double *signal_i,
	  const double *certainty,
	  const double *basis_r,
	  const double *basis_i,
	  const double *applicability,
	  double *result_r,
	  double *result_i,
	  double *cout,
	  const int *signal_dimensions,
	  const int *model_dimensions,
	  const int *result_dimensions,
	  int basis_size,
	  const int *start_indices,
	  const int *stop_indices,
	  const double *roi,
	  const char *cout_func_name,
	  mxArray *cout_data_array,
	  int cout_elements,
	  int is_real)
{
  mxArray *arguments[2];
  mxArray *G_matrix;
  mxArray *h_vector;
  mxArray *x_vector;
  double *G_r;
  double *G_i;
  double *h_r;
  double *h_i;
  double *x_r;
  double *x_i;
  
  mxArray *cout_arguments[5];
  mxArray *G0_matrix;
  mxArray *this_cout_array;
  double *G0_r;
  double *G0_i;
  double *this_cout;
  int cout_num_arguments;
  
  int result_index;
  int model_index;
  int signal_index;
  int displacement;
  
  int resultindex;
  
  int modeldimprod;           /* Product of the first n dimensions. */
  int resultdimprod;          /* Product of the first n dimensions. */
  
  double p_r;
  double p_i;
  double ip_r;
  double ip_i;
  
  int k, k1, k2;
  
  if (is_real)
  {
    G_matrix = mxCreateDoubleMatrix(basis_size, basis_size, mxREAL);
    h_vector = mxCreateDoubleMatrix(basis_size, 1, mxREAL);
    G_r = mxGetPr(G_matrix);
    G_i = NULL;
    h_r = mxGetPr(h_vector);
    h_i = NULL;
  }
  else
  {
    G_matrix = mxCreateDoubleMatrix(basis_size, basis_size, mxCOMPLEX);
    h_vector = mxCreateDoubleMatrix(basis_size, 1, mxCOMPLEX);
    G_r = mxGetPr(G_matrix);
    G_i = mxGetPi(G_matrix);
    h_r = mxGetPr(h_vector);
    h_i = mxGetPi(h_vector);
  }
  
  arguments[0] = G_matrix;
  arguments[1] = h_vector;
  
  /* Initialize the first result index and the model index. */
  result_index = start_indices[0];
  model_index = 0;
  
  /* The centers of the basis functions are supposed to be local
   * origins. Compute the corresponding offset. Note that this is
   * integer division.
   */
  displacement = (model_dimensions[0] - 1) / 2;
  
  /* Compute the product of the first n dimensions of the model and
   * result arrays respectively.
   */
  modeldimprod = model_dimensions[0];
  resultdimprod = result_dimensions[0];
  
  if (cout_func_name)
  {
    init_cout_matrices(basis_r, basis_i, applicability, basis_size,
		       modeldimprod, &cout_num_arguments, cout_arguments,
		       G_matrix, G0_matrix, h_vector, cout_data_array,
		       is_real);
  }
  
  /* Loop over the signal dimensions. */
  for (result_index = start_indices[0];
       result_index <= stop_indices[0];
       result_index++)
  {
    /* If we have a mask region of interest and this is zero, skip
     * this point. (The result array is initialized to zero so we
     * don't need to put anything there ourselves.)
     */
    if (roi && roi[result_index] == 0.0)
      continue;
    
    /* Double loop over the basis functions to compute inner
     * products. Inner products between basis functions and signal
     * are computed when k2 == basis_size.
     */
    for (k1 = 0; k1 < basis_size; k1++)
    {
      for (k2 = k1; k2 <= basis_size; k2++)
      {
	/* Reset the accumulated inner product. */
	ip_r = 0.0;
	ip_i = 0.0;
	
	/* Loop over the dimensions for the basis functions. */
	for (model_index = 0;
	     model_index < model_dimensions[0];
	     model_index++)
	{
	  /* Compute the signal index corresponding to the
	   * current result index and model index.
	   */
	  signal_index = (result_index
			  + model_index
			  - displacement);
	  /* Check if we are outside the signal boundary. It
	   * is implied that the certainty is zero then.
	   */
	  if (signal_index < 0
	      || signal_index >= signal_dimensions[0])
	    continue;
	  
	  if (is_real)
	  {
	    if (k2 == basis_size)
	      p_r = signal_r[signal_index];
	    else
	      p_r = basis_r[model_index + modeldimprod * k2];
	    
	    p_r *= certainty[signal_index];
	    p_r *= applicability[model_index];
	    p_r *= basis_r[model_index + modeldimprod * k1];
	    ip_r += p_r;
	  }
	  else
	  {
	    double ca, b_r, b_i;
	    if (k2 == basis_size)
	    {
	      p_r = signal_r[signal_index];
	      p_i = signal_i[signal_index];
	    }
	    else
	    {
	      p_r = basis_r[model_index + modeldimprod * k2];
	      p_i = basis_i[model_index + modeldimprod * k2];
	    }
	    
	    ca = (certainty[signal_index] * applicability[model_index]);
	    p_r *= ca;
	    p_i *= ca;
	    b_r = basis_r[model_index + modeldimprod * k1];
	    b_i = basis_i[model_index + modeldimprod * k1];
	    /* Notice that this is multiplication with the
	     * conjugate of b.
	     */
	    ip_r += p_r * b_r + p_i * b_i;
	    ip_i += p_i * b_r - p_r * b_i;
	  }
	}

	if (is_real)
	{
	  if (k2 == basis_size)
	    h_r[k1] = ip_r;
	  else
	  {
	    G_r[k1 + k2 * basis_size] = ip_r;
	    G_r[k2 + k1 * basis_size] = ip_r;
	  }
	}
	else
	{
	  if (k2 == basis_size)
	  {
	    h_r[k1] = ip_r;
	    h_i[k1] = ip_i;
	  }
	  else
	  {
	    G_r[k1 + k2 * basis_size] = ip_r;
	    G_i[k1 + k2 * basis_size] = ip_i;
	    G_r[k2 + k1 * basis_size] = ip_r;
	    G_i[k2 + k1 * basis_size] = -ip_i;
	  }
	}
      }
    }
    
    resultindex = (result_index - start_indices[0]);
    
    if (basis_size == 1 && !cout_func_name)
    {
      if (is_real)
	result_r[resultindex] = h_r[0] / G_r[0];
      else
      {
	double w = 1.0 / (G_r[0] * G_r[0] + G_i[0] * G_i[0]);
	result_r[resultindex] = w * (h_r[0] * G_r[0] + h_i[0] * G_i[0]);
	result_i[resultindex] = w * (h_i[0] * G_r[0] - h_r[0] * G_i[0]);
      }
    }
    else
    {
      mexCallMATLAB(1, &x_vector, 2, arguments, "\\");
      if (is_real)
      {
	x_r = mxGetPr(x_vector);
	for (k = 0; k < basis_size; k++)
	  result_r[resultindex + resultdimprod * k] = x_r[k];
      }
      else
      {
	x_r = mxGetPr(x_vector);
	x_i = mxGetPi(x_vector);
	for (k = 0; k < basis_size; k++)
	{
	  result_r[resultindex + resultdimprod * k] = x_r[k];
	  if (x_i != NULL)
	    result_i[resultindex + resultdimprod * k] = x_i[k];
	  else /* Not sure whether this can happen. */
	    result_i[resultindex + resultdimprod * k] = 0.0;
	}
      }
    }
    
    if (cout_func_name)
    {
      cout_arguments[3] = x_vector;
      mexCallMATLAB(1, &this_cout_array, cout_num_arguments, cout_arguments,
		    cout_func_name);
      this_cout = mxGetPr(this_cout_array);
      
      for (k = 0; k < cout_elements; k++)
	cout[resultindex + resultdimprod * k] = this_cout[k];
      
      mxDestroyArray(this_cout_array);
    }
    
    if (basis_size > 1 || cout_func_name)
      mxDestroyArray(x_vector);
  }
}

/* Two signal dimensions, multiple basis functions. */
static void
normconv2(const double *signal_r,
	  const double *signal_i,
	  const double *certainty,
	  const double *basis_r,
	  const double *basis_i,
	  const double *applicability,
	  double *result_r,
	  double *result_i,
	  double *cout,
	  const int *signal_dimensions,
	  const int *model_dimensions,
	  const int *result_dimensions,
	  int basis_size,
	  const int *start_indices,
	  const int *stop_indices,
	  const double *roi,
	  const char *cout_func_name,
	  mxArray *cout_data_array,
	  int cout_elements,
	  int is_real)
{
  mxArray *arguments[2];
  mxArray *G_matrix;
  mxArray *h_vector;
  mxArray *x_vector;
  double *G_r;
  double *G_i;
  double *h_r;
  double *h_i;
  double *x_r;
  double *x_i;
  
  mxArray *cout_arguments[5];
  mxArray *G0_matrix;
  mxArray *this_cout_array;
  double *G0_r;
  double *G0_i;
  double *this_cout;
  int cout_num_arguments;
  
  int result_indices[2];
  int model_indices[2];
  int signal_indices[2];
  int displacements[2];
  
  int signalindex;
  int modelindex;
  int resultindex;
  
  int modeldimprod;           /* Product of the first n dimensions. */
  int resultdimprod;          /* Product of the first n dimensions. */
  
  double p_r;
  double p_i;
  double ip_r;
  double ip_i;
  
  int i, k, k1, k2;
  
  if (is_real)
  {
    G_matrix = mxCreateDoubleMatrix(basis_size, basis_size, mxREAL);
    h_vector = mxCreateDoubleMatrix(basis_size, 1, mxREAL);
    G_r = mxGetPr(G_matrix);
    G_i = NULL;
    h_r = mxGetPr(h_vector);
    h_i = NULL;
  }
  else
  {
    G_matrix = mxCreateDoubleMatrix(basis_size, basis_size, mxCOMPLEX);
    h_vector = mxCreateDoubleMatrix(basis_size, 1, mxCOMPLEX);
    G_r = mxGetPr(G_matrix);
    G_i = mxGetPi(G_matrix);
    h_r = mxGetPr(h_vector);
    h_i = mxGetPi(h_vector);
  }
  
  arguments[0] = G_matrix;
  arguments[1] = h_vector;
  
  /* Initialize the first 2 result indices and the model indices. */
  for (i = 0; i < 2; i++)
  {
    result_indices[i] = start_indices[i];
    model_indices[i] = 0;
  }
  
  /* The centers of the basis functions are supposed to be local
   * origins. Compute the corresponding offsets.
   */
  for (i = 0; i < 2; i++)
  {
    /* Note that this is integer division. */
    displacements[i] = (model_dimensions[i] - 1) / 2;
  }
  
  /* Compute the product of the first n dimensions of the model and
   * result arrays respectively.
   */
  modeldimprod = 1;
  resultdimprod = 1;
  for (i = 0; i < 2; i++)
  {
    modeldimprod *= model_dimensions[i];
    resultdimprod *= result_dimensions[i];
  }
  
  if (cout_func_name)
  {
    init_cout_matrices(basis_r, basis_i, applicability, basis_size,
		       modeldimprod, &cout_num_arguments, cout_arguments,
		       G_matrix, G0_matrix, h_vector, cout_data_array,
		       is_real);
  }
  
  /* Loop over the signal dimensions. */
  for (result_indices[1] = start_indices[1];
       result_indices[1] <= stop_indices[1];
       result_indices[1]++)
  {
    for (result_indices[0] = start_indices[0];
	 result_indices[0] <= stop_indices[0];
	 result_indices[0]++)
    {
      /* If we have a mask region of interest and this is zero, skip
       * this point. (The result array is initialized to zero so we
       * don't need to put anything there ourselves.)
       */
      if (roi)
      {
	int roi_index = (result_indices[1] * signal_dimensions[0]
			 + result_indices[0]);
	if (roi[roi_index] == 0.0)
	  continue;
      }
      
      /* Double loop over the basis functions to compute inner
       * products. Inner products between basis functions and
       * signal are computed when k2 == basis_size.
       */
      for (k1 = 0; k1 < basis_size; k1++)
      {
	for (k2 = k1; k2 <= basis_size; k2++)
	{
	  /* Reset the accumulated inner product. */
	  ip_r = 0.0;
	  ip_i = 0.0;
	  
	  /* Loop over the dimensions for the basis functions. */
	  for (model_indices[1] = 0;
	       model_indices[1] < model_dimensions[1];
	       model_indices[1]++)
	  {
	    /* Compute the signal index corresponding to
	     * the current result index and model index.
	     */
	    signal_indices[1] = (result_indices[1]
				 + model_indices[1]
				 - displacements[1]);
	    /* Check if we are outside the signal
	     * boundary. It is implied that the certainty
	     * is zero then.
	     */
	    if (signal_indices[1] < 0
		|| signal_indices[1] >= signal_dimensions[1])
	      continue;
	    
	    /* Loop over the dimensions for the basis functions. */
	    for (model_indices[0] = 0;
		 model_indices[0] < model_dimensions[0];
		 model_indices[0]++)
	    {
	      signal_indices[0] = (result_indices[0]
				   + model_indices[0]
				   - displacements[0]);
	      if (signal_indices[0] < 0
		  || signal_indices[0] >= signal_dimensions[0])
		continue;
	      
	      signalindex = (signal_indices[1]
			     * signal_dimensions[0]
			     + signal_indices[0]);
	      modelindex = (model_indices[1]
			    * model_dimensions[0]
			    + model_indices[0]);
	      
	      if (is_real)
	      {
		if (k2 == basis_size)
		  p_r = signal_r[signalindex];
		else
		  p_r = basis_r[modelindex + modeldimprod * k2];
		
		p_r *= certainty[signalindex];
		p_r *= applicability[modelindex];
		p_r *= basis_r[modelindex + modeldimprod * k1];
		ip_r += p_r;
	      }
	      else
	      {
		double ca, b_r, b_i;
		if (k2 == basis_size)
		{
		  p_r = signal_r[signalindex];
		  p_i = signal_i[signalindex];
		}
		else
		{
		  p_r = basis_r[modelindex + modeldimprod * k2];
		  p_i = basis_i[modelindex + modeldimprod * k2];
		}
		
		ca = (certainty[signalindex]
		      * applicability[modelindex]);
		p_r *= ca;
		p_i *= ca;
		b_r = basis_r[modelindex + modeldimprod * k1];
		b_i = basis_i[modelindex + modeldimprod * k1];
		/* Notice that this is multiplication with the
		 * conjugate of b.
		 */
		ip_r += p_r * b_r + p_i * b_i;
		ip_i += p_i * b_r - p_r * b_i;
	      }
	    }
	  }
	  
	  if (is_real)
	  {
	    if (k2 == basis_size)
	      h_r[k1] = ip_r;
	    else
	    {
	      G_r[k1 + k2 * basis_size] = ip_r;
	      G_r[k2 + k1 * basis_size] = ip_r;
	    }
	  }
	  else
	  {
	    if (k2 == basis_size)
	    {
	      h_r[k1] = ip_r;
	      h_i[k1] = ip_i;
	    }
	    else
	    {
	      G_r[k1 + k2 * basis_size] = ip_r;
	      G_i[k1 + k2 * basis_size] = ip_i;
	      G_r[k2 + k1 * basis_size] = ip_r;
	      G_i[k2 + k1 * basis_size] = -ip_i;
	    }
	  }
	}
      }
      
      resultindex = ((result_indices[1] - start_indices[1])
		     * result_dimensions[0]
		     + (result_indices[0] - start_indices[0]));
      
      if (basis_size == 1 && !cout_func_name)
      {
	if (is_real)
	  result_r[resultindex] = h_r[0] / G_r[0];
	else
	{
	  double w = 1./(G_r[0] * G_r[0] + G_i[0] * G_i[0]);
	  result_r[resultindex] = w * (h_r[0] * G_r[0] + h_i[0] * G_i[0]);
	  result_i[resultindex] = w * (h_i[0] * G_r[0] - h_r[0] * G_i[0]);
	}
      }
      else
      {
	mexCallMATLAB(1, &x_vector, 2, arguments, "\\");
	if (is_real)
	{
	  x_r = mxGetPr(x_vector);
	  for (k = 0; k < basis_size; k++)
	    result_r[resultindex + resultdimprod * k] = x_r[k];
	}
	else
	{
	  x_r = mxGetPr(x_vector);
	  x_i = mxGetPi(x_vector);
	  for (k = 0; k < basis_size; k++)
	  {
	    result_r[resultindex + resultdimprod * k] = x_r[k];
	    if (x_i != NULL)
	      result_i[resultindex + resultdimprod * k] = x_i[k];
	    else /* Not sure whether this can happen. */
	      result_i[resultindex + resultdimprod * k] = 0.0;
	  }
	}
      }
      
      if (cout_func_name)
      {
	cout_arguments[3] = x_vector;
	mexCallMATLAB(1, &this_cout_array, cout_num_arguments, cout_arguments,
		      cout_func_name);
	this_cout = mxGetPr(this_cout_array);
	
	for (k = 0; k < cout_elements; k++)
	  cout[resultindex + resultdimprod * k] = this_cout[k];
	
	mxDestroyArray(this_cout_array);
      }
      
      if (basis_size > 1 || cout_func_name)
	mxDestroyArray(x_vector);
    }
  }
}

/* Three signal dimensions, multiple basis functions. */
static void
normconv3(const double *signal_r,
	  const double *signal_i,
	  const double *certainty,
	  const double *basis_r,
	  const double *basis_i,
	  const double *applicability,
	  double *result_r,
	  double *result_i,
	  double *cout,
	  const int *signal_dimensions,
	  const int *model_dimensions,
	  const int *result_dimensions,
	  int basis_size,
	  const int *start_indices,
	  const int *stop_indices,
	  const double *roi,
	  const char *cout_func_name,
	  mxArray *cout_data_array,
	  int cout_elements,
	  int is_real)
{
  mxArray *arguments[2];
  mxArray *G_matrix;
  mxArray *h_vector;
  mxArray *x_vector;
  double *G_r;
  double *G_i;
  double *h_r;
  double *h_i;
  double *x_r;
  double *x_i;
  
  mxArray *cout_arguments[5];
  mxArray *G0_matrix;
  mxArray *this_cout_array;
  double *G0_r;
  double *G0_i;
  double *this_cout;
  int cout_num_arguments;
  
  int result_indices[3];
  int model_indices[3];
  int signal_indices[3];
  int displacements[3];
  
  int signalindex;
  int modelindex;
  int resultindex;
  
  int modeldimprod;           /* Product of the first n dimensions. */
  int resultdimprod;          /* Product of the first n dimensions. */
  
  double p_r;
  double p_i;
  double ip_r;
  double ip_i;
  
  int i, k, k1, k2;
  
  if (is_real)
  {
    G_matrix = mxCreateDoubleMatrix(basis_size, basis_size, mxREAL);
    h_vector = mxCreateDoubleMatrix(basis_size, 1, mxREAL);
    G_r = mxGetPr(G_matrix);
    G_i = NULL;
    h_r = mxGetPr(h_vector);
    h_i = NULL;
  }
  else
  {
    G_matrix = mxCreateDoubleMatrix(basis_size, basis_size, mxCOMPLEX);
    h_vector = mxCreateDoubleMatrix(basis_size, 1, mxCOMPLEX);
    G_r = mxGetPr(G_matrix);
    G_i = mxGetPi(G_matrix);
    h_r = mxGetPr(h_vector);
    h_i = mxGetPi(h_vector);
  }
  
  arguments[0] = G_matrix;
  arguments[1] = h_vector;
  
  /* Initialize the first 3 result indices and the model indices. */
  for (i = 0; i < 3; i++)
  {
    result_indices[i] = start_indices[i];
    model_indices[i] = 0;
  }
  
  /* The centers of the basis functions are supposed to be local
   * origins. Compute the corresponding offsets.
   */
  for (i = 0; i < 3; i++)
  {
    /* Note that this is integer division. */
    displacements[i] = (model_dimensions[i] - 1) / 2;
  }
  
  /* Compute the product of the first n dimensions of the model and
   * result arrays respectively.
   */
  modeldimprod = 1;
  resultdimprod = 1;
  for (i = 0; i < 3; i++)
  {
    modeldimprod *= model_dimensions[i];
    resultdimprod *= result_dimensions[i];
  }
  
  if (cout_func_name)
  {
    init_cout_matrices(basis_r, basis_i, applicability, basis_size,
		       modeldimprod, &cout_num_arguments, cout_arguments,
		       G_matrix, G0_matrix, h_vector, cout_data_array,
		       is_real);
  }
  
  /* Loop over the signal dimensions */
  for (result_indices[2] = start_indices[2];
       result_indices[2] <= stop_indices[2];
       result_indices[2]++)
  {
    for (result_indices[1] = start_indices[1];
	 result_indices[1] <= stop_indices[1];
	 result_indices[1]++)
    {
      for (result_indices[0] = start_indices[0];
	   result_indices[0] <= stop_indices[0];
	   result_indices[0]++)
      {
	/* If we have a mask region of interest and this is
	 * zero, skip this point. (The result array is
	 * initialized to zero so we don't need to put
	 * anything there ourselves.)
	 */
	if (roi)
	{
	  int roi_index = ((result_indices[2]
			    * signal_dimensions[1]
			    + result_indices[1])
			   * signal_dimensions[0]
			   + result_indices[0]);
	  if (roi[roi_index] == 0.0)
	    continue;
	}
	
	/* Double loop over the basis functions to compute
	 * inner products. Inner products between basis
	 * functions and signal are computed when k2 ==
	 * basis_size.
	 */
	for (k1 = 0; k1 < basis_size; k1++)
	{
	  for (k2 = k1; k2 <= basis_size; k2++)
	  {
	    /* Reset the accumulated inner product. */
	    ip_r = 0.0;
	    ip_i = 0.0;
	    
	    /* Loop over the dimensions for the basis functions. */
	    for (model_indices[2] = 0;
		 model_indices[2] < model_dimensions[2];
		 model_indices[2]++)
	    {
	      /* Compute the signal index corresponding
	       * to the current result index and model
	       * index.
	       */
	      signal_indices[2] = (result_indices[2]
				   + model_indices[2]
				   - displacements[2]);
	      /* Check if we are outside the signal
	       * boundary. It is implied that the
	       * certainty is zero then.
	       */
	      if (signal_indices[2] < 0
		  || signal_indices[2] >= signal_dimensions[2])
		continue;
	      
	      for (model_indices[1] = 0;
		   model_indices[1] < model_dimensions[1];
		   model_indices[1]++)
	      {
		signal_indices[1] = (result_indices[1]
				     + model_indices[1]
				     - displacements[1]);
		if (signal_indices[1] < 0
		    || (signal_indices[1]
			>= signal_dimensions[1]))
		  continue;
		
		for (model_indices[0] = 0;
		     model_indices[0] < model_dimensions[0];
		     model_indices[0]++)
		{
		  signal_indices[0] = (result_indices[0]
				       + model_indices[0]
				       - displacements[0]);
		  if (signal_indices[0] < 0
		      || (signal_indices[0]
			  >= signal_dimensions[0]))
		    continue;
		  
		  signalindex = ((signal_indices[2]
				  * signal_dimensions[1]
				  + signal_indices[1])
				 * signal_dimensions[0]
				 + signal_indices[0]);
		  modelindex = ((model_indices[2]
				 * model_dimensions[1]
				 + model_indices[1])
				* model_dimensions[0]
				+ model_indices[0]);
		  
		  if (is_real)
		  {
		    if (k2 == basis_size)
		      p_r = signal_r[signalindex];
		    else
		      p_r = basis_r[modelindex + modeldimprod * k2];
		    
		    p_r *= certainty[signalindex];
		    p_r *= applicability[modelindex];
		    p_r *= basis_r[modelindex + modeldimprod * k1];
		    ip_r += p_r;
		  }
		  else
		  {
		    double ca, b_r, b_i;
		    if (k2 == basis_size)
		    {
		      p_r = signal_r[signalindex];
		      p_i = signal_i[signalindex];
		    }
		    else
		    {
		      p_r = basis_r[modelindex + modeldimprod * k2];
		      p_i = basis_i[modelindex + modeldimprod * k2];
		    }
		    
		    ca = (certainty[signalindex]
			  * applicability[modelindex]);
		    p_r *= ca;
		    p_i *= ca;
		    b_r = basis_r[modelindex + modeldimprod * k1];
		    b_i = basis_i[modelindex + modeldimprod * k1];
		    /* Notice that this is multiplication with the
		     * conjugate of b.
		     */
		    ip_r += p_r * b_r + p_i * b_i;
		    ip_i += p_i * b_r - p_r * b_i;
		  }
		}
	      }
	    }
	    
	    if (is_real)
	    {
	      if (k2 == basis_size)
		h_r[k1] = ip_r;
	      else
	      {
		G_r[k1 + k2 * basis_size] = ip_r;
		G_r[k2 + k1 * basis_size] = ip_r;
	      }
	    }
	    else
	    {
	      if (k2 == basis_size)
	      {
		h_r[k1] = ip_r;
		h_i[k1] = ip_i;
	      }
	      else
	      {
		G_r[k1 + k2 * basis_size] = ip_r;
		G_i[k1 + k2 * basis_size] = ip_i;
		G_r[k2 + k1 * basis_size] = ip_r;
		G_i[k2 + k1 * basis_size] = -ip_i;
	      }
	    }
	  }
	}
	
	resultindex = (((result_indices[2] - start_indices[2])
			* result_dimensions[1]
			+ (result_indices[1] - start_indices[1]))
		       * result_dimensions[0]
		       + (result_indices[0] - start_indices[0]));
	
	if (basis_size == 1 && !cout_func_name)
	{
	  if (is_real)
	    result_r[resultindex] = h_r[0] / G_r[0];
	  else
	  {
	    double w = 1./(G_r[0] * G_r[0] + G_i[0] * G_i[0]);
	    result_r[resultindex] = w * (h_r[0] * G_r[0] + h_i[0] * G_i[0]);
	    result_i[resultindex] = w * (h_i[0] * G_r[0] - h_r[0] * G_i[0]);
	  }
	}
	else
	{
	  mexCallMATLAB(1, &x_vector, 2, arguments, "\\");
	  if (is_real)
	  {
	    x_r = mxGetPr(x_vector);
	    for (k = 0; k < basis_size; k++)
	      result_r[resultindex + resultdimprod * k] = x_r[k];
	  }
	  else
	  {
	    x_r = mxGetPr(x_vector);
	    x_i = mxGetPi(x_vector);
	    for (k = 0; k < basis_size; k++)
	    {
	      result_r[resultindex + resultdimprod * k] = x_r[k];
	      if (x_i != NULL)
		result_i[resultindex + resultdimprod * k] = x_i[k];
	      else /* Not sure whether this can happen. */
		result_i[resultindex + resultdimprod * k] = 0.0;
	    }
	  }
	}
	
	if (cout_func_name)
	{
	  cout_arguments[3] = x_vector;
	  mexCallMATLAB(1, &this_cout_array, cout_num_arguments,
			cout_arguments, cout_func_name);
	  this_cout = mxGetPr(this_cout_array);
	  
	  for (k = 0; k < cout_elements; k++)
	    cout[resultindex + resultdimprod * k] = this_cout[k];
	  
	  mxDestroyArray(this_cout_array);
	}
	
	if (basis_size > 1 || cout_func_name)
	  mxDestroyArray(x_vector);
      }
    }
  }
}


static int
calc_single_subscript(const int *dimensions, const int *indices, int dim)
{
  int i;
  int sub = 0;
  for (i = dim - 1; i >= 0; i--)
  {
    sub += indices[i];
    if (i > 0)
      sub *= dimensions[i-1];
  }
  return sub;
}

static int
calc_single_subscript2(const int *dimensions, const int *indices,
		       const int *start_indices, int dim)
{
  int i;
  int sub = 0;
  for (i = dim - 1; i >= 0; i--)
  {
    sub += (indices[i] - start_indices[i]);
    if (i > 0)
      sub *= dimensions[i-1];
  }
  return sub;
}

/* General case, no optimization. */
static void
normconv(const double *signal_r,
	 const double *signal_i,
	 const double *certainty,
	 const double *basis_r,
	 const double *basis_i,
	 const double *applicability,
	 double *result_r,
	 double *result_i,
	 double *cout,
	 int dimensionality,
	 const int *signal_dimensions,
	 const int *model_dimensions,
	 const int *result_dimensions,
	 int basis_size,
	 const int *start_indices,
	 const int *stop_indices,
	 const double *roi,
	 const char *cout_func_name,
	 mxArray *cout_data_array,
	 int cout_elements,
	 int is_real)
{
  mxArray *arguments[2];
  mxArray *G_matrix;
  mxArray *h_vector;
  mxArray *x_vector;
  double *G_r;
  double *G_i;
  double *h_r;
  double *h_i;
  double *x_r;
  double *x_i;
  
  mxArray *cout_arguments[5];
  mxArray *G0_matrix;
  mxArray *this_cout_array;
  double *G0_r;
  double *G0_i;
  double *this_cout;
  int cout_num_arguments;
  
  int *result_indices;
  int *model_indices;
  int *signal_indices;
  int *displacements;
  
  int signalindex;
  int modelindex;
  int resultindex;
  
  int modeldimprod;           /* Product of the first n dimensions. */
  int resultdimprod;          /* Product of the first n dimensions. */
  
  double p_r;
  double p_i;
  double ip_r;
  double ip_i;
  
  int i, k, k1, k2;
  
  int outside;
  
  if (is_real)
  {
    G_matrix = mxCreateDoubleMatrix(basis_size, basis_size, mxREAL);
    h_vector = mxCreateDoubleMatrix(basis_size, 1, mxREAL);
    G_r = mxGetPr(G_matrix);
    G_i = NULL;
    h_r = mxGetPr(h_vector);
    h_i = NULL;
  }
  else
  {
    G_matrix = mxCreateDoubleMatrix(basis_size, basis_size, mxCOMPLEX);
    h_vector = mxCreateDoubleMatrix(basis_size, 1, mxCOMPLEX);
    G_r = mxGetPr(G_matrix);
    G_i = mxGetPi(G_matrix);
    h_r = mxGetPr(h_vector);
    h_i = mxGetPi(h_vector);
  }

  arguments[0] = G_matrix;
  arguments[1] = h_vector;
  
  result_indices = (int *)mxCalloc(dimensionality, sizeof(int));
  model_indices  = (int *)mxCalloc(dimensionality, sizeof(int));
  signal_indices = (int *)mxCalloc(dimensionality, sizeof(int));
  displacements  = (int *)mxCalloc(dimensionality, sizeof(int));
  
  /* Initialize the first n result indices. Note that the model
   * indices already are initialized to zero by the call to
   * mxCalloc.
   */
  for (i = 0; i < dimensionality; i++)
    result_indices[i] = start_indices[i];
  
  /* The centers of the basis functions are supposed to be local
   * origins. Compute the corresponding offsets.
   */
  for (i = 0; i < dimensionality; i++)
  {
    /* Note that this is integer division. */
    displacements[i] = (model_dimensions[i] - 1) / 2;
  }
  
  /* Compute the product of the first n dimensions of the model and
   * result arrays respectively.
   */
  modeldimprod = 1;
  resultdimprod = 1;
  for (i = 0; i < dimensionality; i++)
  {
    modeldimprod *= model_dimensions[i];
    resultdimprod *= result_dimensions[i];
  }
  
  if (cout_func_name)
  {
    init_cout_matrices(basis_r, basis_i, applicability, basis_size,
		       modeldimprod, &cout_num_arguments, cout_arguments,
		       G_matrix, G0_matrix, h_vector, cout_data_array,
		       is_real);
  }
  
  /* Loop over the signal dimensions */
  while (1)
  {
    /* (The indices are incremented at the end of the loop.) */
    
#if DEBUG
    for (i = 0; i < dimensionality; i++)
      mexPrintf("%3d ", result_indices[i]);
    mexPrintf("\n");
#endif
    
    /* If we have a mask region of interest and this is
     * zero, skip this point. (The result array is
     * initialized to zero so we don't need to put
     * anything there ourselves.)
     */
    if (roi)
    {
      int roi_index = calc_single_subscript(signal_dimensions,
					    result_indices,
					    dimensionality);
      if (roi[roi_index] == 0.0)
	goto next_point;
    }
    
    /* Double loop over the basis functions to compute inner
     * products. Inner products between basis functions and signal
     * are computed when k2 == basis_size.
     */
    for (k1 = 0; k1 < basis_size; k1++)
    {
      for (k2 = k1; k2 <= basis_size; k2++)
      {
	/* Reset the accumulated inner product. */		
	ip_r = 0.0;
	ip_i = 0.0;
	
	/* Loop over the dimensions for the basis functions. */
	while (1)
	{
#if DEBUG
	  mexPrintf("Model:");
	  for (i = 0; i < dimensionality; i++)
	    mexPrintf("%3d ", model_indices[i]);
	  mexPrintf("\n");
#endif
	  
	  /* Compute the signal indices corresponding to the
	   * current result indices and model indices.
	   */
	  for (i = 0; i < dimensionality; i++)
	  {
	    signal_indices[i] = result_indices[i]+
	      model_indices[i] - displacements[i];
	  }
	  
	  /* Check if we are outside the signal boundary. It
	   * is implied that the certainty is zero then.
	   */
	  outside = 0;
	  for (i = 0; i < dimensionality; i++)
	  {
	    if (signal_indices[i] < 0
		|| signal_indices[i] >= signal_dimensions[i])
	    {
	      outside = 1;
	      break;
	    }
	  }
	  if (!outside)
	  {
	    signalindex = calc_single_subscript(signal_dimensions,
						signal_indices,
						dimensionality);
	    modelindex = calc_single_subscript(model_dimensions,
					       model_indices,
					       dimensionality);
	    
	    if (is_real)
	    {
	      if (k2 == basis_size)
		p_r = signal_r[signalindex];
	      else
		p_r = basis_r[modelindex + modeldimprod * k2];
	      
	      p_r *= certainty[signalindex];
	      p_r *= applicability[modelindex];
	      p_r *= basis_r[modelindex + modeldimprod * k1];
	      ip_r += p_r;
	    }
	    else
	    {
	      double ca, b_r, b_i;
	      if (k2 == basis_size)
	      {
		p_r = signal_r[signalindex];
		p_i = signal_i[signalindex];
	      }
	      else
	      {
		p_r = basis_r[modelindex + modeldimprod * k2];
		p_i = basis_i[modelindex + modeldimprod * k2];
	      }
	      
	      ca = (certainty[signalindex]
		    * applicability[modelindex]);
	      p_r *= ca;
	      p_i *= ca;
	      b_r = basis_r[modelindex + modeldimprod * k1];
	      b_i = basis_i[modelindex + modeldimprod * k1];
	      /* Notice that this is multiplication with the
	       * conjugate of b.
	       */
	      ip_r += p_r * b_r + p_i * b_i;
	      ip_i += p_i * b_r - p_r * b_i;
	    }
	  }
	  
	  /* Increment the indices. */
	  for (i = 0; i < dimensionality; i++)
	  {
	    model_indices[i]++;
	    if (model_indices[i] >= model_dimensions[i])
	      model_indices[i] = 0;
	    else
	      break;
	  }
	  if (i == dimensionality)
	  {
	    /* Loop finished. By the way, the indices have
	     * all been reset to zero.
	     */
	    break;
	  }
	}
	
	if (is_real)
	{
	  if (k2 == basis_size)
	    h_r[k1] = ip_r;
	  else
	  {
	    G_r[k1 + k2 * basis_size] = ip_r;
	    G_r[k2 + k1 * basis_size] = ip_r;
	  }
	}
	else
	{
	  if (k2 == basis_size)
	  {
	    h_r[k1] = ip_r;
	    h_i[k1] = ip_i;
	  }
	  else
	  {
	    G_r[k1 + k2 * basis_size] = ip_r;
	    G_i[k1 + k2 * basis_size] = ip_i;
	    G_r[k2 + k1 * basis_size] = ip_r;
	    G_i[k2 + k1 * basis_size] = -ip_i;
	  }
	}
      }
    }
    
    resultindex = calc_single_subscript2(result_dimensions,
					 result_indices,
					 start_indices,
					 dimensionality);
	
    mexCallMATLAB(1, &x_vector, 2, arguments, "\\");
    if (is_real)
    {
      x_r = mxGetPr(x_vector);
      for (k = 0; k < basis_size; k++)
	result_r[resultindex + resultdimprod * k] = x_r[k];
    }
    else
    {
      x_r = mxGetPr(x_vector);
      x_i = mxGetPi(x_vector);
      for (k = 0; k < basis_size; k++)
      {
	result_r[resultindex + resultdimprod * k] = x_r[k];
	if (x_i != NULL)
	  result_i[resultindex + resultdimprod * k] = x_i[k];
	else /* Not sure whether this can happen. */
	  result_i[resultindex + resultdimprod * k] = 0.0;
      }
    }
    
    if (cout_func_name)
    {
      cout_arguments[3] = x_vector;
      mexCallMATLAB(1, &this_cout_array,
		    cout_num_arguments, cout_arguments, cout_func_name);
      this_cout = mxGetPr(this_cout_array);
      
      for (k = 0; k < cout_elements; k++)
	cout[resultindex + resultdimprod * k] = this_cout[k];
      
      mxDestroyArray(this_cout_array);
    }
    
    mxDestroyArray(x_vector);
    
    /* Increment the indices. */
    next_point:
    for (i = 0; i < dimensionality; i++)
    {
      result_indices[i]++;
      if (result_indices[i] > stop_indices[i])
	result_indices[i] = start_indices[i];
      else
	break;
    }
    if (i == dimensionality)
    {
      /* Loop finished */	
      break;
    }
  }
  
  mxFree(result_indices);
  mxFree(model_indices);
  mxFree(signal_indices);
  mxFree(displacements);
}

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int i;
  int dimensionality;
  int sn, cn, bn, an;
  const int *sdim = NULL;
  const int *cdim = NULL;
  const int *bdim = NULL;
  const int *adim = NULL;
  int *signal_dimensions = NULL;
  int *model_dimensions = NULL;
  int *result_dimensions = NULL;
  int *cout_dimensions = NULL;
  mxArray *result_array = NULL;
  mxArray *cout_array = NULL;
  double *cout = NULL;
  double *region_of_interest = NULL;
  int startindex, stopindex;
  int *start_indices = NULL;
  int *stop_indices = NULL;
  int basis_size;
  int optimization = 1;
  const mxArray *roi_array = NULL;
  double *roi = NULL;
  const mxArray *options_array = NULL;
  mxArray *cout_data_array = NULL;
  char *cout_func_name = NULL;
  int cout_elements = 0;
  int roi_type = NO_ROI;
  int is_real = 1;
  const double *signal_i;
  const double *basis_i;
  
  /* Check the number of input and output arguments. */
  if (nrhs < 4)
    mexErrMsgTxt("Too few input arguments.");
  if (nrhs > 6)
    mexErrMsgTxt("Too many input arguments.");
  if (nlhs > 2)
    mexErrMsgTxt("Too many output arguments.");
  
  /* Check the formats of the first four input arguments. */
  for (i = 0; i < 4; i++)
  {
    if (!mxIsNumeric(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsDouble(prhs[i]))
      mexErrMsgTxt("The first four input arguments must be full numeric arrays, stored as doubles.");
    
    if ((i == 1 || i == 3) && mxIsComplex(prhs[i]))
      mexErrMsgTxt("Certainty and applicability must be real.");

    if ((i == 0 || i == 2) && mxIsComplex(prhs[i]))
      is_real = 0;

    /* We won't deal with empty arrays. */
    if (mxIsEmpty(prhs[i]))
      mexErrMsgTxt("Empty arrays not allowed.");
  }
  
  /* Look for region_of_interest and/or options. */
  roi_array = NULL;
  options_array = NULL;
  for (i = 4; i < nrhs; i++)
  {
    if (mxIsStruct(prhs[i]))
    {
      if (options_array)
	mexErrMsgTxt("Syntax error, two struct arrays unexpected.");
      options_array = prhs[i];
    }
    else if (mxIsNumeric(prhs[i]) && !mxIsComplex(prhs[i])
	     && !mxIsSparse(prhs[i]) && mxIsDouble(prhs[i]))
    {
      if (roi_array)
	mexErrMsgTxt("Syntax error, two region_of_interest unexpected.");
      roi_array = prhs[i];
    }
    else
      mexErrMsgTxt("Argument looks neither like region_of_interest, nor like options.");
  }
  
  
  /* Get the dimensionalities. Trailing singleton dimensions are
   * ignored.
   */
  
  /* signal */
  sdim = mxGetDimensions(prhs[0]);
  for (sn = mxGetNumberOfDimensions(prhs[0]); sn > 0; sn--)
    if (sdim[sn-1] > 1)
      break;
  
  /* certainty */
  cdim = mxGetDimensions(prhs[1]);
  for (cn = mxGetNumberOfDimensions(prhs[1]); cn > 0; cn--)
    if (cdim[cn-1] > 1)
      break;
  
  /* basis */
  bdim = mxGetDimensions(prhs[2]);
  for (bn = mxGetNumberOfDimensions(prhs[2]); bn > 0; bn--)
    if (bdim[bn-1] > 1)
      break;
  
  /* applicability */
  adim = mxGetDimensions(prhs[3]);
  for (an = mxGetNumberOfDimensions(prhs[3]); an > 0; an--)
    if (adim[an-1] > 1)
      break;
  
  /* Let the dimensionality of the convolution be the largest of the
   * given parameters.
   */
  dimensionality = sn;
  if (cn > dimensionality)
    dimensionality = cn;
  if (bn-1 > dimensionality)
    dimensionality = bn-1;
  if (an > dimensionality)
    dimensionality = an;
  
  /* Fix to manage the rather uninteresting special case of all
   * parameters being scalars.
   */
  if (dimensionality == 0)
    dimensionality = 1;
  
  /* Build dimension vectors for signal and model. Also check the
   * consistency of the dimensions.
   */
  signal_dimensions = (int *)mxCalloc(dimensionality, sizeof(int));
  model_dimensions  = (int *)mxCalloc(dimensionality, sizeof(int));
  
  if (cn != sn)
    mexErrMsgTxt("Signal and certainty must have the same size.");
  for (i = 0; i < sn; i++)
    if (cdim[i] != sdim[i])
      mexErrMsgTxt("Signal and certainty must have the same size.");
  
  for (i = 0; i < sn; i++)
    signal_dimensions[i] = sdim[i];
  for (; i < dimensionality; i++)
    signal_dimensions[i] = 1;
  
  /* The interpretation of the dimensions of basis is made in a very
   * friendly manner. It is required that the non-singleton
   * dimensions of applicability are present at the beginning of
   * basis. If all the following dimensions are singleton it is
   * assumed to be only one basis function. If exactly one of the
   * following dimensions is larger than one, that is assumed to be
   * the number of basis functions. Otherwise an error is declared.
   */
  if (bn < an)
    mexErrMsgTxt("Each basis function must have the same size as applicability.");
  
  for (i = 0; i < an; i++)
    if (bdim[i] != adim[i])
      mexErrMsgTxt("Each basis function must have the same size as applicability.");
  
  for (; i < bn; i++)
    if (bdim[i] > 1)
      break;
  
  if (i == bn)
    basis_size = 1;
  else
  {
    basis_size = bdim[i];
    for (i++; i < bn; i++)
      if (bdim[i] > 1)
	break;
    if (i < bn)
      mexErrMsgTxt("Each basis function must have the same size as applicability.");
  }
  
  for (i = 0; i < an; i++)
    model_dimensions[i] = adim[i];
  for (; i < dimensionality; i++)
    model_dimensions[i] = 1;
  
  
  /* Check the options parameter, if present. */
  if (options_array)
  {
    /* Look for region_of_interest. */
    const mxArray *opt = mxGetField(options_array, 0,
				    "region_of_interest");
    if (opt)
    {
      if (!mxIsNumeric(opt) || mxIsComplex(opt) || !mxIsDouble(opt))
	mexErrMsgTxt("The region_of_interest field in options must be a real and full numeric array, stored as doubles.");
      
      if (roi_array)
	mexErrMsgTxt("Cannot have both a region_of_interest parameter and a region_of_interest field in options.");
      else
	roi_array = opt;
    }
    
    /* Look for non_optimized. */
    opt = mxGetField(options_array, 0, "non_optimized"); 
    if (opt && mxGetScalar(opt) != 0.0)
      optimization = 0;
    
    /* Look for cout_func. */
    opt = mxGetField(options_array, 0, "cout_func");
    if (opt)
    {
      if (!mxIsChar(opt))
	mexErrMsgTxt("The cout_func field of options is expected to be a string.");
      cout_func_name = mxCalloc(mxGetN(opt) + 1, 1);
      mxGetString(opt, cout_func_name, mxGetN(opt) + 1);
    }
    
    /* Look for cout_data. */
    cout_data_array = mxGetField(options_array, 0, "cout_data");
  }
  
  /* Check the validity of the region of interest. */
  if (roi_array)
  {
    if (mxGetNumberOfDimensions(roi_array) == 2
	&& mxGetM(roi_array) == sn
	&& mxGetN(roi_array) == 2)
    {
      region_of_interest = mxGetPr(roi_array);
      for (i = 0; i < sn; i++)
      {
	startindex = region_of_interest[i];
	stopindex = region_of_interest[i + sn];
	if (startindex < 1
	    || startindex > stopindex
	    || stopindex > sdim[i])
	  mexErrMsgTxt("Invalid region of interest.");
      }
      roi_type = BOX_ROI;
    }
    else
    {
      int ok = 1;
      const int *roidim = mxGetDimensions(roi_array);
      int roin;
      
      for (roin = mxGetNumberOfDimensions(roi_array); roin > 0; roin--)
	if (roidim[roin-1] > 1)
	  break;
      
      if (roin != sn)
	ok = 0;

      for (i = 0; i < sn; i++)
	if (roidim[i] != sdim[i])
	  ok = 0;

      if (!ok)
	mexErrMsgTxt("Region of interest must be either an N by 2 matrix, where N is the dimensionality of the signal, or have the same size as the signal.");

      roi_type = MASK_ROI;
      roi = mxGetPr(roi_array);
    }
  }
  
  /* Create the start and stop indices. */
  start_indices = (int *)mxCalloc(dimensionality, sizeof(int));
  stop_indices  = (int *)mxCalloc(dimensionality, sizeof(int));
  
  if (roi_type != BOX_ROI)
  {
    for (i = 0; i < sn; i++)
    {
      start_indices[i] = 0;
      stop_indices[i] = sdim[i] - 1;
    }
  }
  else
  {
    region_of_interest = mxGetPr(roi_array);
    for (i = 0; i < sn; i++)
    {
      start_indices[i] = region_of_interest[i] - 1;
      stop_indices[i] = region_of_interest[i + sn] - 1;
    }
  }
  
  for (i = sn; i < dimensionality; i++)
  {
    start_indices[i] = 0;
    stop_indices[i] = 0;
  }
  
  /* Compute the output dimensions. */
  result_dimensions = (int *)mxCalloc(dimensionality + 1, sizeof(int));
  for (i = 0; i < dimensionality; i++)
    result_dimensions[i] = stop_indices[i] - start_indices[i] + 1;
  result_dimensions[dimensionality] = basis_size;
  
  /* Create the r array. */
  result_array = mxCreateNumericArray(dimensionality+1, result_dimensions,
				      mxDOUBLE_CLASS,
				      is_real ? mxREAL : mxCOMPLEX);
  
  if (nlhs > 1 && !cout_func_name)
    mexErrMsgTxt("Too many output arguments. Cout requires a cout_func field in the options parameter.");
  
  if (nlhs == 1 && cout_func_name)
  {
    /* Nowhere to put the output certainty. Pointless to compute it. */
    cout_func_name = NULL;
  }
  
  /* Create the cout array. */
  if (nlhs == 2)
  {
    mxArray *arguments[5];
    mxArray *G;
    mxArray *h;
    mxArray *sample_cout;
    int n_arguments = 4;
    int cout_num_dims;
    const int *cout_dims;
    /* First call cout_func to see how large output certainty we get.
     * We need one identity matrix of size MxM, where M is the number
     * of basis functions, plus one Mx1 zero vector.
     */
    G = mxCreateDoubleMatrix(basis_size, basis_size, mxREAL);
    for (i = 0; i < basis_size; i++)
      mxGetPr(G)[i + i * basis_size] = 1.0;
    h = mxCreateDoubleMatrix(basis_size, 1, mxREAL);
    arguments[0] = G;
    arguments[1] = G;
    arguments[2] = h;
    arguments[3] = h;
    if (cout_data_array)
    {
      arguments[4] = cout_data_array;
      n_arguments = 5;
    }
    
    mexCallMATLAB(1, &sample_cout, n_arguments, arguments,
		  cout_func_name);
    cout_num_dims = mxGetNumberOfDimensions(sample_cout);
    cout_dims = mxGetDimensions(sample_cout);
    cout_dimensions = (int *)mxCalloc(dimensionality + cout_num_dims,
				      sizeof(int));
    for (i = 0; i < dimensionality; i++)
      cout_dimensions[i] = result_dimensions[i];
    for (; i < dimensionality + cout_num_dims; i++)
      cout_dimensions[i] = cout_dims[i - dimensionality];
    
    cout_elements = mxGetM(sample_cout) * mxGetN(sample_cout);
    
    mxDestroyArray(sample_cout);
    
    cout_array = mxCreateNumericArray(dimensionality + cout_num_dims,
				      cout_dimensions, mxDOUBLE_CLASS, mxREAL);
    cout = mxGetPr(cout_array);
  }
  
  if (!is_real)
  {
    signal_i = mxGetPi(prhs[0]);
    if (!signal_i)
      signal_i = mxCalloc(mxGetM(prhs[0]) * mxGetN(prhs[0]),
			  sizeof(double));
    basis_i = mxGetPi(prhs[2]);
    if (!basis_i)
      basis_i = mxCalloc(mxGetM(prhs[2]) * mxGetN(prhs[2]),
			 sizeof(double));
  }
  
  /* Call the appropriate computational routine. */
  if (dimensionality > 3)
    optimization = 0;
  
  switch (dimensionality * optimization)
  {
    case 0: /* General case, no optimization */
      normconv(mxGetPr(prhs[0]),
	       signal_i,
	       mxGetPr(prhs[1]),
	       mxGetPr(prhs[2]),
	       basis_i,
	       mxGetPr(prhs[3]),
	       mxGetPr(result_array),
	       mxGetPi(result_array),
	       cout,
	       dimensionality,
	       signal_dimensions,
	       model_dimensions,
	       result_dimensions,
	       basis_size,
	       start_indices,
	       stop_indices,
	       roi,
	       cout_func_name,
	       cout_data_array,
	       cout_elements,
	       is_real);
      break;
      
    case 1: /* One dimension. */
      normconv1(mxGetPr(prhs[0]),
		signal_i,
		mxGetPr(prhs[1]),
		mxGetPr(prhs[2]),
		basis_i,
		mxGetPr(prhs[3]),
		mxGetPr(result_array),
		mxGetPi(result_array),
		cout,
		signal_dimensions,
		model_dimensions,
		result_dimensions,
		basis_size,
		start_indices,
		stop_indices,
		roi,
		cout_func_name,
		cout_data_array,
		cout_elements,
		is_real);
      break;
      
    case 2: /* Two dimensions. */
      normconv2(mxGetPr(prhs[0]),
		signal_i,
		mxGetPr(prhs[1]),
		mxGetPr(prhs[2]),
		basis_i,
		mxGetPr(prhs[3]),
		mxGetPr(result_array),
		mxGetPi(result_array),
		cout,
		signal_dimensions,
		model_dimensions,
		result_dimensions,
		basis_size,
		start_indices,
		stop_indices,
		roi,
		cout_func_name,
		cout_data_array,
		cout_elements,
		is_real);
      break;
      
    case 3: /* Three dimensions. */
      normconv3(mxGetPr(prhs[0]),
		signal_i,
		mxGetPr(prhs[1]),
		mxGetPr(prhs[2]),
		basis_i,
		mxGetPr(prhs[3]),
		mxGetPr(result_array),
		mxGetPi(result_array),
		cout,
		signal_dimensions,
		model_dimensions,
		result_dimensions,
		basis_size,
		start_indices,
		stop_indices,
		roi,
		cout_func_name,
		cout_data_array,
		cout_elements,
		is_real);
      break;
      
    default:
      mexErrMsgTxt("Internal error, impossible case.");
      break;
  }
  
  /* Free allocated memory. */
  mxFree(start_indices);
  mxFree(stop_indices);
  
  /* Output the computed result. */
  plhs[0] = result_array;
  if (nrhs > 1)
    plhs[1] = cout_array;
}
