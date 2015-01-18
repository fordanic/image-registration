#include "mex.h"
#include <string.h>

#define RECURSION_SIZE_LIMIT 4

/*
 * ANTIGRADIENT
 * 
 * See the help text in antigradient.m for a description.
 * 
 * Author: Gunnar Farnebäck
 *         Medical Informatics
 *         Linköping University, Sweden
 *         gunnar@imt.liu.se
 */


/*************** 2D ****************/

void
solve_directly2D(double *f, double *rhs, double *f_out,
		 int M, int N)
{
  int s = M * N;
  int dims[2];
  mxArray *A_array;
  mxArray *b_array;
  double *A;
  double *b;
  int k;
  int i, j;
  mxArray *x_array;
  mxArray *input_arrays[2];
  
  dims[0] = s;
  dims[1] = s;
  A_array = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  A = mxGetPr(A_array);
  
  dims[1] = 1;
  b_array = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  b = mxGetPr(b_array);
  
  k = 0;
  for (j = 0; j < N; j++)
    for (i = 0; i < M; i++)
    {
      A[k + s * (j*M+i)] = -4;
      if (i > 0)
	A[k + s * (j*M+i-1)] = 1 + (i == M-1);
      if (i < M-1)
	A[k + s * (j*M+i+1)] = 1 + (i == 0);
      if (j > 0)
	A[k + s * ((j-1)*M+i)] = 1 + (j == N-1);
      if (j < N-1)
	A[k + s * ((j+1)*M+i)] = 1 + (j == 0);
      
      b[k] = rhs[i+j*M];
      
      k++;
    }
  
  for (i = 0; i < s*s; i++)
    A[i] += 1.0 / (s*s);
  
  input_arrays[0] = A_array;
  input_arrays[1] = b_array;
  mexCallMATLAB(1, &x_array, 2, input_arrays, "\\");
  memcpy(f_out, mxGetPr(x_array), s * sizeof(*f));
  mxDestroyArray(x_array);
  mxDestroyArray(A_array);
  mxDestroyArray(b_array);
}


/* Gauss-Seidel smoothing iteration. Red-black ordering. */
void
gauss_seidel2D(double *f, double *d, int M, int N)
{
  int pass;
  int i, j;
  int index;
  
  for (pass = 0; pass <= 1; pass++)
  {
    for (j = 0; j < N; j++)
      for (i = 0; i < M; i++)
      {
	double new_f;
	if ((i + j) % 2 != pass)
	  continue;
	
	index = i + j * M;
	new_f = -d[index];
	if (i == 0)
	  new_f += f[index + 1];
	else
	  new_f += f[index - 1];
	if (i == M-1)
	  new_f += f[index - 1];
	else
	  new_f += f[index + 1];
	
	if (j == 0)
	  new_f += f[index + M];
	else
	  new_f += f[index - M];
	if (j == N-1)
	  new_f += f[index - M];
	else
	  new_f += f[index + M];
	
	f[index] = 0.25 * new_f;
      }
  }
}


void
downsample2D(double *rhs, int M, int N,
	     double *rhs_coarse, int Mhalf, int Nhalf)
{
  int i, j;
  int index1;
  int index2;
  
  if (M % 2 == 0 && N % 2 == 0)
  {
    for (j = 0; j < Nhalf; j++)
      for (i = 0; i < Mhalf; i++)
      {
	index1 = (j * Mhalf + i);
	index2 = (2 * j * M + 2 * i);
	
	rhs_coarse[index1] = (rhs[index2]
			      + rhs[index2 + M]
			      + rhs[index2 + 1]
			      + rhs[index2 + M + 1]);
      }
  }
  
  if (M % 2 == 1 && N % 2 == 0)
  {
    for (j = 0; j < Nhalf; j++)
    {
      index1 = (j * Mhalf);
      index2 = (2 * j * M);
      rhs_coarse[index1] = (rhs[index2]
			    + rhs[index2 + M]
			    + rhs[index2 + 1]
			    + rhs[index2 + M + 1]);
      for (i = 1; i < Mhalf - 1; i++)
      {
	index1++;
	index2 += 2;
	rhs_coarse[index1] = (rhs[index2]
			      + rhs[index2 + M]
			      + 0.5 * (rhs[index2 - 1]
				       + rhs[index2 + 1]
				       + rhs[index2 + M - 1]
				       + rhs[index2 + M + 1]));
      }
      index1++;
      index2 += 2;
      rhs_coarse[index1] = (rhs[index2]
			    + rhs[index2 + M]
			    + rhs[index2 - 1]
			    + rhs[index2 + M - 1]);
    }
  }
  
  if (M % 2 == 0 && N % 2 == 1)
  {
    for (i = 0; i < Mhalf; i++)
    {
      index1 = i;
      index2 = 2 * i;
      rhs_coarse[index1] = (rhs[index2]
			    + rhs[index2 + 1]
			    + rhs[index2 + M]
			    + rhs[index2 + M + 1]);
    }
    
    for (j = 1; j < Nhalf - 1; j++)
    {
      for (i = 0; i < Mhalf; i++)
      {
	index1 = (j * Mhalf + i);
	index2 = (2 * j * M + 2 * i);
	rhs_coarse[index1] = (rhs[index2]
			      + rhs[index2 + 1]
			      + 0.5 * (rhs[index2 - M]
				       + rhs[index2 + M]
				       + rhs[index2 - M + 1]
				       + rhs[index2 + M + 1]));
      }
    }
    
    for (i = 0; i < Mhalf; i++)
    {
      index1 = (j * Mhalf + i);
      index2 = (2 * j * M + 2 * i);
      rhs_coarse[index1] = (rhs[index2]
			    + rhs[index2 + 1]
			    + rhs[index2 - M]
			    + rhs[index2 - M + 1]);
    }
  }
  
  if (M % 2 == 1 && N % 2 == 1)
  {
    for (j = 1; j < Nhalf - 1; j++)
    {
      for (i = 1; i < Mhalf - 1; i++)
      {
	index1 = (j * Mhalf + i);
	index2 = (2 * j * M + 2 * i);
	rhs_coarse[index1] = (rhs[index2]
			      + 0.5 * (rhs[index2 - 1]
				       + rhs[index2 + 1]
				       + rhs[index2 - M]
				       + rhs[index2 + M])
			      + 0.25* (rhs[index2 - M - 1]
				       + rhs[index2 - M + 1]
				       + rhs[index2 + M - 1]
				       + rhs[index2 + M + 1]));
      }
    }
    
    for (i = 1; i < Mhalf - 1; i++)
    {
      index1 = i;
      index2 = 2 * i;
      rhs_coarse[index1] = (rhs[index2]
			    + rhs[index2 + M]
			    + 0.5 * (rhs[index2 - 1]
				     + rhs[index2 + 1]
				     + rhs[index2 + M - 1]
				     + rhs[index2 + M + 1]));
      index1 = ((Nhalf - 1) * Mhalf + i);
      index2 = (2 * (Nhalf - 1) * M + 2 * i);
      rhs_coarse[index1] = (rhs[index2]
			    + rhs[index2 - M]
			    + 0.5 * (rhs[index2 - 1]
				     + rhs[index2 + 1]
				     + rhs[index2 - M - 1]
				     + rhs[index2 - M + 1]));
    }
    
    for (j = 1; j < Nhalf - 1; j++)
    {
      index1 = (j * Mhalf);
      index2 = (2 * j * M);
      rhs_coarse[index1] = (rhs[index2]
			    + rhs[index2 + 1]
			    + 0.5 * (rhs[index2 - M]
				     + rhs[index2 + M]
				     + rhs[index2 - M + 1]
				     + rhs[index2 + M + 1]));
      index1 = (j * Mhalf + Mhalf - 1);
      index2 = (2 * j * M + 2 * (Mhalf - 1));
      rhs_coarse[index1] = (rhs[index2]
			    + rhs[index2 - 1]
			    + 0.5 * (rhs[index2 - M]
				     + rhs[index2 + M]
				     + rhs[index2 - M - 1]
				     + rhs[index2 + M - 1]));
    }
    
    index1 = 0;
    index2 = 0;
    rhs_coarse[index1] = (rhs[index2]
			  + rhs[index2 + 1]
			  + rhs[index2 + M]
			  + rhs[index2 + M + 1]);
    
    index1 = (Nhalf - 1) * Mhalf;
    index2 = 2 * (Nhalf - 1) * M;
    rhs_coarse[index1] = (rhs[index2]
			  + rhs[index2 + 1]
			  + rhs[index2 - M]
			  + rhs[index2 - M + 1]);
    
    index1 = (Mhalf - 1);
    index2 = (2 * (Mhalf - 1));
    rhs_coarse[index1] = (rhs[index2]
			  + rhs[index2 - 1]
			  + rhs[index2 + M]
			  + rhs[index2 + M - 1]);
    
    index1 = ((Nhalf - 1) * Mhalf + Mhalf - 1);
    index2 = (2 * (Nhalf - 1) * M + 2 * (Mhalf - 1));
    rhs_coarse[index1] = (rhs[index2]
			  + rhs[index2 - 1]
			  + rhs[index2 - M]
			  + rhs[index2 - M - 1]);
  }
}


/* Upsample and apply correction. Bilinear interpolation. */
void
upsample2D(double *rhs, int M, int N,
	   double *v, int Mhalf, int Nhalf,
	   double *f_out)
{
  int i, j;
  int index1, index2;
  double ce, no, so, we, ea, nw, sw, ne, se;
  int CE, NO, SO, WE, EA, SW, NE, SE;
  
  if (M % 2 == 0 && N % 2 == 0)
  {
    for (j = 0; j < Nhalf; j++)
      for (i = 0; i < Mhalf; i++)
      {
	index1 = (j * Mhalf + i);
	index2 = (2 * j * M + 2 * i);
	ce = v[index1];
	no = v[index1 - 1];
	so = v[index1 + 1];
	we = v[index1 - Mhalf];
	ea = v[index1 + Mhalf];
	nw = v[index1 - Mhalf - 1];
	sw = v[index1 - Mhalf + 1];
	ne = v[index1 + Mhalf - 1];
	se = v[index1 + Mhalf + 1];
	CE = index2;
	SO = index2 + 1;
	EA = index2 + M;
	SE = index2 + M + 1;
	
	/* Fine pixels northwest of coarse pixel center. */
	if (i == 0)
	{
	  if (j == 0) /* NW corner. */
	    f_out[CE] += ce - 0.25 * rhs[CE];
	  else        /* North edge. */
	    f_out[CE] += 0.75 * ce + 0.25 * (we - rhs[CE]);
	}
	else
	{
	  if (j == 0) /* West edge. */
	    f_out[CE] += 0.75 * ce + 0.25 * (no - rhs[CE]);
	  else        /* Inner point. */
	    f_out[CE] += (0.5625 * ce + 0.1875 * (no + we) + 0.0625 * nw);
	}
	
	/* Fine pixels southwest of coarse pixel center. */
	if (i == Mhalf - 1)
	{
	  if (j == 0) /* SW corner. */
	    f_out[SO] += ce - 0.25 * rhs[SO];
	  else        /* South edge. */
	    f_out[SO] += 0.75 * ce  + 0.25 * (we - rhs[SO]);
	}
	else
	{
	  if (j == 0) /* West edge. */
	    f_out[SO] += 0.75 * ce + 0.25 * (so - rhs[SO]);
	  else        /* Inner point. */
	    f_out[SO] += 0.5625 * ce + 0.1875 * (so + we) + 0.0625 * sw;
	}
	
	/* Fine pixels northeast of coarse pixel center. */
	if (i == 0)
	{
	  if (j == Nhalf - 1) /* NE corner. */
	    f_out[EA] += ce - 0.25 * rhs[EA];
	  else              /* North edge. */
	    f_out[EA] += 0.75 * ce + 0.25 * (ea - rhs[EA]);
	}
	else
	{
	  if (j == Nhalf - 1) /* East edge. */
	    f_out[EA] += 0.75 * ce + 0.25 * (no - rhs[EA]);
	  else              /* Inner point. */
	    f_out[EA] += 0.5625 * ce + 0.1875 * (no + ea) + 0.0625 * ne;
	}
	
	/* Fine pixels southeast of coarse pixel center. */
	if (i == Mhalf - 1)
	{
	  if (j == Nhalf - 1) /* SE corner. */
	    f_out[SE] += ce - 0.25 * rhs[SE];
	  else              /* South edge.*/
	    f_out[SE] += 0.75 * ce + 0.25 * (ea - rhs[SE]);
	}
	else
	{
	  if (j == Nhalf - 1) /* East edge. */
	    f_out[SE] += 0.75 * ce + 0.25 * (so - rhs[SE]);
	  else              /* Inner point. */
	    f_out[SE] += 0.5625 * ce + 0.1875 * (so + ea) + 0.0625 * se;
	}
      }
  }
  
  if (M % 2 == 1 && N % 2 == 0)
  {
    for (j = 0; j < Nhalf; j++)
      for (i = 0; i < Mhalf; i++)
      {
	index1 = (j * Mhalf + i);
	index2 = (2 * j * M + 2 * i);
	ce = v[index1];
	so = v[index1 + 1];
	we = v[index1 - Mhalf];
	ea = v[index1 + Mhalf];
	sw = v[index1 - Mhalf + 1];
	se = v[index1 + Mhalf + 1];
	CE = index2;
	SO = index2 + 1;
	EA = index2 + M;
	SE = index2 + M + 1;
	NO = index2 - 1;
	NE = index2 + M - 1;
	
	/* Fine pixels west of coarse pixel center. */
	if (j == 0) /* West edge. */
	{
	  if (i == 0) /* NW corner */
	    f_out[CE] += ce - 0.125 * (rhs[SO] + rhs[CE] - rhs[EA]);
	  else if (i == Mhalf - 1) /* SW corner */
	    f_out[CE] += ce - 0.125 * (rhs[NO] + rhs[CE] - rhs[EA]);
	  else
	    f_out[CE] += ce - 0.25 * rhs[CE];
	}
	else        /* Inner point. */
	  f_out[CE] += 0.75 * ce + 0.25 * we;
	
	/* Fine pixels southsouthwest of coarse pixel center. */
	if (i < Mhalf - 1)
	{
	  if (j == 0) /* West edge. */
	    f_out[SO] += 0.5 * (ce + so) - 0.25 * rhs[SO];
	  else        /* Inner point. */
	    f_out[SO] += 0.375 * (ce + so) + 0.125 * (we + sw);
	}
	
	/* Fine pixels east of coarse pixel center. */
	if (j == Nhalf - 1) /* East edge. */
	{
	  if (i == 0) /* NE corner */
	    f_out[EA] += ce - 0.125 * (rhs[SE] + rhs[EA] - rhs[CE]);
	  else if (i == Mhalf - 1) /* SE corner */
	    f_out[EA] += ce - 0.125 * (rhs[NE] + rhs[EA] - rhs[CE]);
	  else
	    f_out[EA] += ce - 0.25 * rhs[EA];
	}
	else          /* Inner point. */
	  f_out[EA] += 0.75 * ce + 0.25 * ea;
	
	/* Fine pixels southsoutheast of coarse pixel center. */
	if (i < Mhalf - 1)
	{
	  if (j == Nhalf - 1) /* East edge. */
	    f_out[SE] += 0.5 * (ce + so) - 0.25 * rhs[SE];
	  else                /* Inner point. */
	    f_out[SE] += 0.375 * (ce + so) + 0.125 * (ea + se);
	}
      }
  }
  
  if (M % 2 == 0 && N % 2 == 1)
  {
    for (j = 0; j < Nhalf; j++)
      for (i = 0; i < Mhalf; i++)
      {
	index1 = (j * Mhalf + i);
	index2 = (2 * j * M + 2 * i);
	ce = v[index1];
	so = v[index1 + 1];
	no = v[index1 - 1];
	ea = v[index1 + Mhalf];
	se = v[index1 + Mhalf + 1];
	ne = v[index1 + Mhalf - 1];
	CE = index2;
	SO = index2 + 1;
	EA = index2 + M;
	SE = index2 + M + 1;
	WE = index2 - M;
	SW = index2 - M + 1;
	
	/* Fine pixels north of coarse pixel center. */
	if (i == 0) /* North edge. */
	{
	  if (j == 0) /* NW corner */
	    f_out[CE] += ce - 0.125 * (rhs[EA] + rhs[CE] - rhs[SO]);
	  else if (j == Nhalf - 1) /* NE corner */
	    f_out[CE] += ce - 0.125 * (rhs[WE] + rhs[CE] - rhs[SO]);
	  else
	    f_out[CE] += ce - 0.25 * rhs[CE];
	}
	else        /* Inner point. */
	  f_out[CE] += 0.75 * ce + 0.25 * no;
	
	/* Fine pixels eastnortheast of coarse pixel center. */
	if (j < Nhalf - 1)
	{
	  if (i == 0) /* North edge. */
	    f_out[EA] += 0.5 * (ce + ea) - 0.25 * rhs[EA];
	  else        /* Inner point. */
	    f_out[EA] += 0.375 * (ce + ea) + 0.125 * (no + ne);
	}
	
	/* Fine pixels south of coarse pixel center. */
	if (i == Mhalf - 1) /* South edge. */
	{
	  if (j == 0) /* SW corner */
	    f_out[SO] += ce - 0.125 * (rhs[SE] + rhs[SO] - rhs[CE]);
	  else if (j == Nhalf - 1) /* SE corner */
	    f_out[SO] += ce - 0.125 * (rhs[SW] + rhs[SO] - rhs[CE]);
	  else
	    f_out[SO] += ce - 0.25 * rhs[SO];
	}
	else          /* Inner point. */
	  f_out[SO] += 0.75 * ce + 0.25 * so;
	
	/* Fine pixels eastsoutheast of coarse pixel center. */
	if (j < Nhalf - 1)
	{
	  if (i == Mhalf - 1) /* South edge. */
	    f_out[SE] += 0.5 * (ce + ea) - 0.25 * rhs[SE];
	  else                /* Inner point. */
	    f_out[SE] += 0.375 * (ce + ea) + 0.125 * (so + se);
	}
      }
  }
  
  if (M % 2 == 1 && N % 2 == 1)
  {
    for (j = 0; j < Nhalf; j++)
      for (i = 0; i < Mhalf; i++)
      {
	index1 = (j * Mhalf + i);
	index2 = (2 * j * M + 2 * i);
	ce = v[index1];
	so = v[index1 + 1];
	ea = v[index1 + Mhalf];
	se = v[index1 + Mhalf + 1];
	CE = index2;
	SO = index2 + 1;
	EA = index2 + M;
	SE = index2 + M + 1;
	
	/* Fine pixels on coarse pixel center. */
	f_out[CE] += ce;
	
	/* Fine pixels south of coarse pixel center. */
	if (i < Mhalf - 1)
	  f_out[SO] += 0.5 * (ce + so);
	
	/* Fine pixels east of coarse pixel center. */
	if (j < Nhalf - 1) /* South edge. */
	  f_out[EA] += 0.5 * (ce + ea);
	
	/* Fine pixels southeast of coarse pixel center. */
	if (i < Mhalf - 1 && j < Nhalf - 1)
	  f_out[SE] += 0.25 * (ce + ea + so + se);
      }
  }
}


/* Recursive multigrid function.*/
void
poisson_multigrid2D(double *f, double *d,
		    int n1, int n2, int nm,
		    double *f_out,
		    int M, int N, int *directly_solved)
{
  int i, j;
  int k;
  double *r;
  double *r_downsampled;
  double *v;
  int Mhalf;
  int Nhalf;
  
  /* Solve a sufficiently small problem directly. */
  if (M < RECURSION_SIZE_LIMIT || N < RECURSION_SIZE_LIMIT)
  {
    solve_directly2D(f, d, f_out, M, N);
    *directly_solved = 1;
    return;
  }
  *directly_solved = 0;
  
  /* Initialize solution. */
  memcpy(f_out, f, M * N * sizeof(*f_out));
  
  /* Pre-smoothing. */
  for (k = 0; k < n1; k++)
    gauss_seidel2D(f_out, d, M, N);
  
  /* Compute residual. */
  r = mxCalloc(M * N, sizeof(*r));
  for (j = 0; j < N; j++)
    for (i = 0; i < M; i++)
    {
      int index = j * M + i;
      double residual = d[index] + 4 * f_out[index];
      if (i == 0)
	residual -= f_out[index + 1];
      else
	residual -= f_out[index - 1];
      
      if (i == M - 1)
	residual -= f_out[index - 1];
      else
	residual -= f_out[index + 1];
      
      if (j == 0)
	residual -= f_out[index + M];
      else
	residual -= f_out[index - M];
      
      if (j == N - 1)
	residual -= f_out[index - M];
      else
	residual -= f_out[index + M];
      
      r[index] = residual;
    }
  
  /* Downsample residual. */
  Mhalf = (M + 1) / 2;
  Nhalf = (N + 1) / 2;
  r_downsampled = mxCalloc(Mhalf * Nhalf, sizeof(*r_downsampled));
  downsample2D(r, M, N, r_downsampled, Mhalf, Nhalf);
  
  /* Recurse to compute a correction. */
  v = mxCalloc(Mhalf * Nhalf, sizeof(*v));
  for (k = 0; k < nm; k++)
  {
    int directly_solved;
    poisson_multigrid2D(v, r_downsampled, n1, n2, nm, v, Mhalf, Nhalf,
			&directly_solved);
    if (directly_solved)
      break;
  }
  
  upsample2D(r, M, N, v, Mhalf, Nhalf, f_out);
  
  /* Post-smoothing. */
  for (k = 0; k < n2; k++)
    gauss_seidel2D(f_out, d, M, N);
  
  mxFree(r);
  mxFree(r_downsampled);
  mxFree(v);
}


/* It is assumed that f_out is initialized to zero when called. */
void
poisson_full_multigrid2D(double *rhs, int number_of_iterations,
			 int M, int N, double *f_out)
{
  double *rhs_downsampled;
  double *f_coarse;
  int k;
  
  /* Unless already coarsest scale, first recurse to coarser scale. */
  if (M >= RECURSION_SIZE_LIMIT && N >= RECURSION_SIZE_LIMIT)
  {
    /* Downsample right hand side. */
    int Mhalf = (M + 1) / 2;
    int Nhalf = (N + 1) / 2;
    rhs_downsampled = mxCalloc(Mhalf * Nhalf, sizeof(*rhs_downsampled));
    downsample2D(rhs, M, N, rhs_downsampled, Mhalf, Nhalf);
    
    f_coarse = mxCalloc(Mhalf * Nhalf, sizeof(*f_coarse));
    poisson_full_multigrid2D(rhs_downsampled, number_of_iterations,
			     Mhalf, Nhalf, f_coarse);
    
    /* Upsample the coarse result. */
    upsample2D(rhs, M, N, f_coarse, Mhalf, Nhalf, f_out);

    mxFree(f_coarse);
    mxFree(rhs_downsampled);
  }
  
  /* Perform number_of_iterations standard multigrid cycles. */
  for (k = 0; k < number_of_iterations; k++)
  {
    int directly_solved;
    poisson_multigrid2D(f_out, rhs, 2, 2, 2, f_out, M, N,
			&directly_solved);
    if (directly_solved)
      break;
  }
}


void
antigradient2D(double *g, double mu, int number_of_iterations,
	       int M, int N, double *f_out)
{
  double *rhs;
  double sum;
  double mean;
  int i, j;
  
  /* Compute right hand side of Poisson problem with Neumann
   * boundary conditions, discretized by finite differences.
   */
  rhs = mxCalloc(M * N, sizeof(*rhs));
  for (j = 0; j < N; j++)
    for (i = 0; i < M; i++)
    {
      int index1 = j * M + i;
      int index2 = index1 + M * N;
      double d = 0.0;
      
      if (i == 0)
	d = g[index1 + 1] + g[index1];
      else if (i == M - 1)
	d = - g[index1] - g[index1 - 1];
      else
	d = 0.5 * (g[index1 + 1] - g[index1 - 1]);
      
      if (j == 0)
	d += g[index2 + M] + g[index2];
      else if (j == N - 1)
	d += - g[index2] - g[index2 - M];
      else
	d += 0.5 * (g[index2 + M] - g[index2 - M]);
      rhs[index1] = d;
    }
  
  /* Solve the equation system with the full multigrid algorithm.
   * Use W cycles and 2 presmoothing and 2 postsmoothing
   * Gauss-Seidel iterations.
   */
  poisson_full_multigrid2D(rhs, number_of_iterations, M, N, f_out);
  
  /* Fix the mean value. */
  sum = 0.0;
  for (i = 0; i < M * N; i++)
    sum += f_out[i];
  
  mean = sum / (M * N);
  for (i = 0; i < M * N; i++)
  {
    f_out[i] -= mean;
    f_out[i] += mu;
  }
}


/*************** 3D ****************/

void
solve_directly3D(double *f, double *rhs,
		 double *f_out,
		 int M, int N, int P)
{
  int s = M * N * P;
  int dims[2];
  mxArray *A_array;
  mxArray *b_array;
  double *A;
  double *b;
  int k;
  int i, j, p;
  mxArray *x_array;
  mxArray *input_arrays[2];
  
  dims[0] = s;
  dims[1] = s;
  A_array = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  A = mxGetPr(A_array);
  
  dims[1] = 1;
  b_array = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  b = mxGetPr(b_array);
  
  k = 0;
  for (p = 0; p < P; p++)
    for (j = 0; j < N; j++)
      for (i = 0; i < M; i++)
      {
	A[k + s * ((p*N+j)*M+i)] = -6;
	if (i > 0)
	  A[k + s * ((p*N+j)*M+i-1)] = 1 + (i == M-1);
	if (i < M-1)
	  A[k + s * ((p*N+j)*M+i+1)] = 1 + (i == 0);
	if (j > 0)
	  A[k + s * ((p*N+(j-1))*M+i)] = 1 + (j == N-1);
	if (j < N-1)
	  A[k + s * ((p*N+(j+1))*M+i)] = 1 + (j == 0);
	if (p > 0)
	  A[k + s * (((p-1)*N+j)*M+i)] = 1 + (p == P-1);
	if (p < P-1)
	  A[k + s * (((p+1)*N+j)*M+i)] = 1 + (p == 0);
	
	b[k] = rhs[(p*N+j)*M+i];
	
	k++;
      }
  
  for (i = 0; i < s*s; i++)
    A[i] += 1.0 / (s*s);
  
  input_arrays[0] = A_array;
  input_arrays[1] = b_array;
  mexCallMATLAB(1, &x_array, 2, input_arrays, "\\");
  memcpy(f_out, mxGetPr(x_array), s * sizeof(*f));
  mxDestroyArray(x_array);
  mxDestroyArray(A_array);
  mxDestroyArray(b_array);
}


/* Gauss-Seidel smoothing iteration. Red-black ordering. */
void
gauss_seidel3D(double *f, double *d, int M, int N, int P)
{
  int pass;
  int i, j, p;
  int index;
  int MN = M * N;
  
  for (pass = 0; pass <= 1; pass++)
  {
    for (p = 0; p < P; p++)
      for (j = 0; j < N; j++)
	for (i = 0; i < M; i++)
	{
	  double new_f;
	  if ((i + j + p) % 2 != pass)
	    continue;
	  
	  index = (p * N + j) * M + i;
	  new_f = -d[index];
	  if (i == 0)
	    new_f += f[index + 1];
	  else
	    new_f += f[index - 1];
	  if (i == M-1)
	    new_f += f[index - 1];
	  else
	    new_f += f[index + 1];
	  
	  if (j == 0)
	    new_f += f[index + M];
	  else
	    new_f += f[index - M];
	  if (j == N-1)
	    new_f += f[index - M];
	  else
	    new_f += f[index + M];
	  
	  if (p == 0)
	    new_f += f[index + MN];
	  else
	    new_f += f[index - MN];
	  if (p == P-1)
	    new_f += f[index - MN];
	  else
	    new_f += f[index + MN];
	  
	  f[index] = (1 / 6.0) * new_f;
	}
  }
}

void
downsample3D(double *rhs, int M, int N, int P,
             double *rhs_coarse, int Mhalf, int Nhalf, int Phalf)
{
  int i, j, p;
  int index1;
  int index2;
  int MN = M * N;
  
  if (M % 2 == 0 && N % 2 == 0 && P % 2 == 0)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          rhs_coarse[index1] = (0.5 * (rhs[index2]
                                       + rhs[index2 + 1]
                                       + rhs[index2 + M]
                                       + rhs[index2 + M + 1]
                                       + rhs[index2 + MN]
                                       + rhs[index2 + MN + 1]
                                       + rhs[index2 + MN + M]
                                       + rhs[index2 + MN + M + 1]));
        }
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 0 && P % 2 == 0)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        index1 = (p * Nhalf + j) * Mhalf + 0;
        index2 = (2 * p * N + 2 * j) * M + 2 * 0;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + 1]
                                     + rhs[index2 + M]
                                     + rhs[index2 + M + 1]
                                     + rhs[index2 + MN]
                                     + rhs[index2 + MN + 1]
                                     + rhs[index2 + MN + M]
                                     + rhs[index2 + MN + M + 1]));
        
        for (i = 1; i < Mhalf - 1; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          rhs_coarse[index1] = (0.5 * (rhs[index2]
                                       + rhs[index2 + M]
                                       + rhs[index2 + MN]
                                       + rhs[index2 + MN + M])
                                + 0.25 * (rhs[index2 + 1]
                                          + rhs[index2 - 1]
                                          + rhs[index2 + M + 1]
                                          + rhs[index2 + M - 1]
                                          + rhs[index2 + MN + 1]
                                          + rhs[index2 + MN - 1]
                                          + rhs[index2 + MN + M + 1]
                                          + rhs[index2 + MN + M - 1]));
        }
        
        index1 = (p * Nhalf + j) * Mhalf + (Mhalf - 1);
        index2 = (2 * p * N + 2 * j) * M + 2 * (Mhalf - 1);
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 - 1]
                                     + rhs[index2 + M]
                                     + rhs[index2 + M - 1]
                                     + rhs[index2 + MN]
                                     + rhs[index2 + MN - 1]
                                     + rhs[index2 + MN + M]
                                     + rhs[index2 + MN + M - 1]));
      }
    }
  }
  
  if (M % 2 == 0 && N % 2 == 1 && P % 2 == 0)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (i = 0; i < Mhalf; i++)
      {
        index1 = (p * Nhalf + 0) * Mhalf + i;
        index2 = (2 * p * N + 2 * 0) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + 1]
                                     + rhs[index2 + M]
                                     + rhs[index2 + M + 1]
                                     + rhs[index2 + MN]
                                     + rhs[index2 + MN + 1]
                                     + rhs[index2 + MN + M]
                                     + rhs[index2 + MN + M + 1]));
      }
      
      for (j = 1; j < Nhalf - 1; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          rhs_coarse[index1] = (0.5 * (rhs[index2]
                                       + rhs[index2 + 1]
                                       + rhs[index2 + MN]
                                       + rhs[index2 + MN + 1])
                                + 0.25 * (rhs[index2 + M]
                                          + rhs[index2 + M + 1]
                                          + rhs[index2 - M]
                                          + rhs[index2 - M + 1]
                                          + rhs[index2 + MN + M]
                                          + rhs[index2 + MN + M + 1]
                                          + rhs[index2 + MN - M]
                                          + rhs[index2 + MN - M + 1]));
        }
      }
      
      for (i = 0; i < Mhalf; i++)
      {
        index1 = (p * Nhalf + (Nhalf - 1)) * Mhalf + i;
        index2 = (2 * p * N + 2 * (Nhalf - 1)) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + 1]
                                     + rhs[index2 - M]
                                     + rhs[index2 - M + 1]
                                     + rhs[index2 + MN]
                                     + rhs[index2 + MN + 1]
                                     + rhs[index2 + MN - M]
                                     + rhs[index2 + MN - M + 1]));
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 1 && P % 2 == 0)
  {
    for (p = 0; p < Phalf; p++)
    {
      index1 = (p * Nhalf + 0) * Mhalf + 0;
      index2 = (2 * p * N + 2 * 0) * M + 2 * 0;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 + 1]
                                   + rhs[index2 + M]
                                   + rhs[index2 + M + 1]
                                   + rhs[index2 + MN]
                                   + rhs[index2 + MN + 1]
                                   + rhs[index2 + MN + M]
                                   + rhs[index2 + MN + M + 1]));
      
      for (i = 1; i < Mhalf - 1; i++)
      {
        index1 = (p * Nhalf + 0) * Mhalf + i;
        index2 = (2 * p * N + 2 * 0) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + M]
                                     + rhs[index2 + MN]
                                     + rhs[index2 + MN + M])
                              + 0.25 * (rhs[index2 + 1]
                                        + rhs[index2 - 1]
                                        + rhs[index2 + M + 1]
                                        + rhs[index2 + M - 1]
                                        + rhs[index2 + MN + 1]
                                        + rhs[index2 + MN - 1]
                                        + rhs[index2 + MN + M + 1]
                                        + rhs[index2 + MN + M - 1]));
      }
      
      index1 = (p * Nhalf + 0) * Mhalf + (Mhalf - 1);
      index2 = (2 * p * N + 2 * 0) * M + 2 * (Mhalf - 1);
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 - 1]
                                   + rhs[index2 + M]
                                   + rhs[index2 + M - 1]
                                   + rhs[index2 + MN]
                                   + rhs[index2 + MN - 1]
                                   + rhs[index2 + MN + M]
                                   + rhs[index2 + MN + M - 1]));
      
      for (j = 1; j < Nhalf - 1; j++)
      {
        index1 = (p * Nhalf + j) * Mhalf + 0;
        index2 = (2 * p * N + 2 * j) * M + 2 * 0;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + 1]
                                     + rhs[index2 + MN]
                                     + rhs[index2 + MN + 1])
                              + 0.25 * (rhs[index2 + M]
                                        + rhs[index2 + M + 1]
                                        + rhs[index2 - M]
                                        + rhs[index2 - M + 1]
                                        + rhs[index2 + MN + M]
                                        + rhs[index2 + MN + M + 1]
                                        + rhs[index2 + MN - M]
                                        + rhs[index2 + MN - M + 1]));
        
        for (i = 1; i < Mhalf - 1; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          rhs_coarse[index1] = (0.5 * (rhs[index2]
                                       + rhs[index2 + MN])
                                + 0.25 * (rhs[index2 + 1]
                                          + rhs[index2 - 1]
                                          + rhs[index2 + M]
                                          + rhs[index2 - M]
                                          + rhs[index2 + MN + 1]
                                          + rhs[index2 + MN - 1]
                                          + rhs[index2 + MN + M]
                                          + rhs[index2 + MN - M])
                                + 0.125 * (rhs[index2 + M + 1]
                                           + rhs[index2 + M - 1]
                                           + rhs[index2 - M + 1]
                                           + rhs[index2 - M - 1]
                                           + rhs[index2 + MN + M + 1]
                                           + rhs[index2 + MN + M - 1]
                                           + rhs[index2 + MN - M + 1]
                                           + rhs[index2 + MN - M - 1]));
        }
        
        index1 = (p * Nhalf + j) * Mhalf + (Mhalf - 1);
        index2 = (2 * p * N + 2 * j) * M + 2 * (Mhalf - 1);
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 - 1]
                                     + rhs[index2 + MN]
                                     + rhs[index2 + MN - 1])
                              + 0.25 * (rhs[index2 + M]
                                        + rhs[index2 + M - 1]
                                        + rhs[index2 - M]
                                        + rhs[index2 - M - 1]
                                        + rhs[index2 + MN + M]
                                        + rhs[index2 + MN + M - 1]
                                        + rhs[index2 + MN - M]
                                        + rhs[index2 + MN - M - 1]));
      }
      
      index1 = (p * Nhalf + (Nhalf - 1)) * Mhalf + 0;
      index2 = (2 * p * N + 2 * (Nhalf - 1)) * M + 2 * 0;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 + 1]
                                   + rhs[index2 - M]
                                   + rhs[index2 - M + 1]
                                   + rhs[index2 + MN]
                                   + rhs[index2 + MN + 1]
                                   + rhs[index2 + MN - M]
                                   + rhs[index2 + MN - M + 1]));
      
      for (i = 1; i < Mhalf - 1; i++)
      {
        index1 = (p * Nhalf + (Nhalf - 1)) * Mhalf + i;
        index2 = (2 * p * N + 2 * (Nhalf - 1)) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 - M]
                                     + rhs[index2 + MN]
                                     + rhs[index2 + MN - M])
                              + 0.25 * (rhs[index2 + 1]
                                        + rhs[index2 - 1]
                                        + rhs[index2 - M + 1]
                                        + rhs[index2 - M - 1]
                                        + rhs[index2 + MN + 1]
                                        + rhs[index2 + MN - 1]
                                        + rhs[index2 + MN - M + 1]
                                        + rhs[index2 + MN - M - 1]));
      }
      
      index1 = (p * Nhalf + (Nhalf - 1)) * Mhalf + (Mhalf - 1);
      index2 = (2 * p * N + 2 * (Nhalf - 1)) * M + 2 * (Mhalf - 1);
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 - 1]
                                   + rhs[index2 - M]
                                   + rhs[index2 - M - 1]
                                   + rhs[index2 + MN]
                                   + rhs[index2 + MN - 1]
                                   + rhs[index2 + MN - M]
                                   + rhs[index2 + MN - M - 1]));
    }
  }
  
  if (M % 2 == 0 && N % 2 == 0 && P % 2 == 1)
  {
    for (j = 0; j < Nhalf; j++)
    {
      for (i = 0; i < Mhalf; i++)
      {
        index1 = (0 * Nhalf + j) * Mhalf + i;
        index2 = (2 * 0 * N + 2 * j) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + 1]
                                     + rhs[index2 + M]
                                     + rhs[index2 + M + 1]
                                     + rhs[index2 + MN]
                                     + rhs[index2 + MN + 1]
                                     + rhs[index2 + MN + M]
                                     + rhs[index2 + MN + M + 1]));
      }
    }
    
    for (p = 1; p < Phalf - 1; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          rhs_coarse[index1] = (0.5 * (rhs[index2]
                                       + rhs[index2 + 1]
                                       + rhs[index2 + M]
                                       + rhs[index2 + M + 1])
                                + 0.25 * (rhs[index2 + MN]
                                          + rhs[index2 + MN + 1]
                                          + rhs[index2 + MN + M]
                                          + rhs[index2 + MN + M + 1]
                                          + rhs[index2 - MN]
                                          + rhs[index2 - MN + 1]
                                          + rhs[index2 - MN + M]
                                          + rhs[index2 - MN + M + 1]));
        }
      }
    }
    
    for (j = 0; j < Nhalf; j++)
    {
      for (i = 0; i < Mhalf; i++)
      {
        index1 = ((Phalf - 1) * Nhalf + j) * Mhalf + i;
        index2 = (2 * (Phalf - 1) * N + 2 * j) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + 1]
                                     + rhs[index2 + M]
                                     + rhs[index2 + M + 1]
                                     + rhs[index2 - MN]
                                     + rhs[index2 - MN + 1]
                                     + rhs[index2 - MN + M]
                                     + rhs[index2 - MN + M + 1]));
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 0 && P % 2 == 1)
  {
    for (j = 0; j < Nhalf; j++)
    {
      index1 = (0 * Nhalf + j) * Mhalf + 0;
      index2 = (2 * 0 * N + 2 * j) * M + 2 * 0;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 + 1]
                                   + rhs[index2 + M]
                                   + rhs[index2 + M + 1]
                                   + rhs[index2 + MN]
                                   + rhs[index2 + MN + 1]
                                   + rhs[index2 + MN + M]
                                   + rhs[index2 + MN + M + 1]));
      
      for (i = 1; i < Mhalf - 1; i++)
      {
        index1 = (0 * Nhalf + j) * Mhalf + i;
        index2 = (2 * 0 * N + 2 * j) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + M]
                                     + rhs[index2 + MN]
                                     + rhs[index2 + MN + M])
                              + 0.25 * (rhs[index2 + 1]
                                        + rhs[index2 - 1]
                                        + rhs[index2 + M + 1]
                                        + rhs[index2 + M - 1]
                                        + rhs[index2 + MN + 1]
                                        + rhs[index2 + MN - 1]
                                        + rhs[index2 + MN + M + 1]
                                        + rhs[index2 + MN + M - 1]));
      }
      
      index1 = (0 * Nhalf + j) * Mhalf + (Mhalf - 1);
      index2 = (2 * 0 * N + 2 * j) * M + 2 * (Mhalf - 1);
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 - 1]
                                   + rhs[index2 + M]
                                   + rhs[index2 + M - 1]
                                   + rhs[index2 + MN]
                                   + rhs[index2 + MN - 1]
                                   + rhs[index2 + MN + M]
                                   + rhs[index2 + MN + M - 1]));
    }
    
    for (p = 1; p < Phalf - 1; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        index1 = (p * Nhalf + j) * Mhalf + 0;
        index2 = (2 * p * N + 2 * j) * M + 2 * 0;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + 1]
                                     + rhs[index2 + M]
                                     + rhs[index2 + M + 1])
                              + 0.25 * (rhs[index2 + MN]
                                        + rhs[index2 + MN + 1]
                                        + rhs[index2 + MN + M]
                                        + rhs[index2 + MN + M + 1]
                                        + rhs[index2 - MN]
                                        + rhs[index2 - MN + 1]
                                        + rhs[index2 - MN + M]
                                        + rhs[index2 - MN + M + 1]));
        
        for (i = 1; i < Mhalf - 1; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          rhs_coarse[index1] = (0.5 * (rhs[index2]
                                       + rhs[index2 + M])
                                + 0.25 * (rhs[index2 + 1]
                                          + rhs[index2 - 1]
                                          + rhs[index2 + M + 1]
                                          + rhs[index2 + M - 1]
                                          + rhs[index2 + MN]
                                          + rhs[index2 + MN + M]
                                          + rhs[index2 - MN]
                                          + rhs[index2 - MN + M])
                                + 0.125 * (rhs[index2 + MN + 1]
                                           + rhs[index2 + MN - 1]
                                           + rhs[index2 + MN + M + 1]
                                           + rhs[index2 + MN + M - 1]
                                           + rhs[index2 - MN + 1]
                                           + rhs[index2 - MN - 1]
                                           + rhs[index2 - MN + M + 1]
                                           + rhs[index2 - MN + M - 1]));
        }
        
        index1 = (p * Nhalf + j) * Mhalf + (Mhalf - 1);
        index2 = (2 * p * N + 2 * j) * M + 2 * (Mhalf - 1);
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 - 1]
                                     + rhs[index2 + M]
                                     + rhs[index2 + M - 1])
                              + 0.25 * (rhs[index2 + MN]
                                        + rhs[index2 + MN - 1]
                                        + rhs[index2 + MN + M]
                                        + rhs[index2 + MN + M - 1]
                                        + rhs[index2 - MN]
                                        + rhs[index2 - MN - 1]
                                        + rhs[index2 - MN + M]
                                        + rhs[index2 - MN + M - 1]));
      }
    }
    
    for (j = 0; j < Nhalf; j++)
    {
      index1 = ((Phalf - 1) * Nhalf + j) * Mhalf + 0;
      index2 = (2 * (Phalf - 1) * N + 2 * j) * M + 2 * 0;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 + 1]
                                   + rhs[index2 + M]
                                   + rhs[index2 + M + 1]
                                   + rhs[index2 - MN]
                                   + rhs[index2 - MN + 1]
                                   + rhs[index2 - MN + M]
                                   + rhs[index2 - MN + M + 1]));
      
      for (i = 1; i < Mhalf - 1; i++)
      {
        index1 = ((Phalf - 1) * Nhalf + j) * Mhalf + i;
        index2 = (2 * (Phalf - 1) * N + 2 * j) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + M]
                                     + rhs[index2 - MN]
                                     + rhs[index2 - MN + M])
                              + 0.25 * (rhs[index2 + 1]
                                        + rhs[index2 - 1]
                                        + rhs[index2 + M + 1]
                                        + rhs[index2 + M - 1]
                                        + rhs[index2 - MN + 1]
                                        + rhs[index2 - MN - 1]
                                        + rhs[index2 - MN + M + 1]
                                        + rhs[index2 - MN + M - 1]));
      }
      
      index1 = ((Phalf - 1) * Nhalf + j) * Mhalf + (Mhalf - 1);
      index2 = (2 * (Phalf - 1) * N + 2 * j) * M + 2 * (Mhalf - 1);
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 - 1]
                                   + rhs[index2 + M]
                                   + rhs[index2 + M - 1]
                                   + rhs[index2 - MN]
                                   + rhs[index2 - MN - 1]
                                   + rhs[index2 - MN + M]
                                   + rhs[index2 - MN + M - 1]));
    }
  }
  
  if (M % 2 == 0 && N % 2 == 1 && P % 2 == 1)
  {
    for (i = 0; i < Mhalf; i++)
    {
      index1 = (0 * Nhalf + 0) * Mhalf + i;
      index2 = (2 * 0 * N + 2 * 0) * M + 2 * i;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 + 1]
                                   + rhs[index2 + M]
                                   + rhs[index2 + M + 1]
                                   + rhs[index2 + MN]
                                   + rhs[index2 + MN + 1]
                                   + rhs[index2 + MN + M]
                                   + rhs[index2 + MN + M + 1]));
    }
    
    for (j = 1; j < Nhalf - 1; j++)
    {
      for (i = 0; i < Mhalf; i++)
      {
        index1 = (0 * Nhalf + j) * Mhalf + i;
        index2 = (2 * 0 * N + 2 * j) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + 1]
                                     + rhs[index2 + MN]
                                     + rhs[index2 + MN + 1])
                              + 0.25 * (rhs[index2 + M]
                                        + rhs[index2 + M + 1]
                                        + rhs[index2 - M]
                                        + rhs[index2 - M + 1]
                                        + rhs[index2 + MN + M]
                                        + rhs[index2 + MN + M + 1]
                                        + rhs[index2 + MN - M]
                                        + rhs[index2 + MN - M + 1]));
      }
    }
    
    for (i = 0; i < Mhalf; i++)
    {
      index1 = (0 * Nhalf + (Nhalf - 1)) * Mhalf + i;
      index2 = (2 * 0 * N + 2 * (Nhalf - 1)) * M + 2 * i;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 + 1]
                                   + rhs[index2 - M]
                                   + rhs[index2 - M + 1]
                                   + rhs[index2 + MN]
                                   + rhs[index2 + MN + 1]
                                   + rhs[index2 + MN - M]
                                   + rhs[index2 + MN - M + 1]));
    }
    
    for (p = 1; p < Phalf - 1; p++)
    {
      for (i = 0; i < Mhalf; i++)
      {
        index1 = (p * Nhalf + 0) * Mhalf + i;
        index2 = (2 * p * N + 2 * 0) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + 1]
                                     + rhs[index2 + M]
                                     + rhs[index2 + M + 1])
                              + 0.25 * (rhs[index2 + MN]
                                        + rhs[index2 + MN + 1]
                                        + rhs[index2 + MN + M]
                                        + rhs[index2 + MN + M + 1]
                                        + rhs[index2 - MN]
                                        + rhs[index2 - MN + 1]
                                        + rhs[index2 - MN + M]
                                        + rhs[index2 - MN + M + 1]));
      }
      
      for (j = 1; j < Nhalf - 1; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          rhs_coarse[index1] = (0.5 * (rhs[index2]
                                       + rhs[index2 + 1])
                                + 0.25 * (rhs[index2 + M]
                                          + rhs[index2 + M + 1]
                                          + rhs[index2 - M]
                                          + rhs[index2 - M + 1]
                                          + rhs[index2 + MN]
                                          + rhs[index2 + MN + 1]
                                          + rhs[index2 - MN]
                                          + rhs[index2 - MN + 1])
                                + 0.125 * (rhs[index2 + MN + M]
                                           + rhs[index2 + MN + M + 1]
                                           + rhs[index2 + MN - M]
                                           + rhs[index2 + MN - M + 1]
                                           + rhs[index2 - MN + M]
                                           + rhs[index2 - MN + M + 1]
                                           + rhs[index2 - MN - M]
                                           + rhs[index2 - MN - M + 1]));
        }
      }
      
      for (i = 0; i < Mhalf; i++)
      {
        index1 = (p * Nhalf + (Nhalf - 1)) * Mhalf + i;
        index2 = (2 * p * N + 2 * (Nhalf - 1)) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + 1]
                                     + rhs[index2 - M]
                                     + rhs[index2 - M + 1])
                              + 0.25 * (rhs[index2 + MN]
                                        + rhs[index2 + MN + 1]
                                        + rhs[index2 + MN - M]
                                        + rhs[index2 + MN - M + 1]
                                        + rhs[index2 - MN]
                                        + rhs[index2 - MN + 1]
                                        + rhs[index2 - MN - M]
                                        + rhs[index2 - MN - M + 1]));
      }
    }
    
    for (i = 0; i < Mhalf; i++)
    {
      index1 = ((Phalf - 1) * Nhalf + 0) * Mhalf + i;
      index2 = (2 * (Phalf - 1) * N + 2 * 0) * M + 2 * i;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 + 1]
                                   + rhs[index2 + M]
                                   + rhs[index2 + M + 1]
                                   + rhs[index2 - MN]
                                   + rhs[index2 - MN + 1]
                                   + rhs[index2 - MN + M]
                                   + rhs[index2 - MN + M + 1]));
    }
    
    for (j = 1; j < Nhalf - 1; j++)
    {
      for (i = 0; i < Mhalf; i++)
      {
        index1 = ((Phalf - 1) * Nhalf + j) * Mhalf + i;
        index2 = (2 * (Phalf - 1) * N + 2 * j) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + 1]
                                     + rhs[index2 - MN]
                                     + rhs[index2 - MN + 1])
                              + 0.25 * (rhs[index2 + M]
                                        + rhs[index2 + M + 1]
                                        + rhs[index2 - M]
                                        + rhs[index2 - M + 1]
                                        + rhs[index2 - MN + M]
                                        + rhs[index2 - MN + M + 1]
                                        + rhs[index2 - MN - M]
                                        + rhs[index2 - MN - M + 1]));
      }
    }
    
    for (i = 0; i < Mhalf; i++)
    {
      index1 = ((Phalf - 1) * Nhalf + (Nhalf - 1)) * Mhalf + i;
      index2 = (2 * (Phalf - 1) * N + 2 * (Nhalf - 1)) * M + 2 * i;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 + 1]
                                   + rhs[index2 - M]
                                   + rhs[index2 - M + 1]
                                   + rhs[index2 - MN]
                                   + rhs[index2 - MN + 1]
                                   + rhs[index2 - MN - M]
                                   + rhs[index2 - MN - M + 1]));
    }
  }
  
  if (M % 2 == 1 && N % 2 == 1 && P % 2 == 1)
  {
    index1 = (0 * Nhalf + 0) * Mhalf + 0;
    index2 = (2 * 0 * N + 2 * 0) * M + 2 * 0;
    rhs_coarse[index1] = (0.5 * (rhs[index2]
                                 + rhs[index2 + 1]
                                 + rhs[index2 + M]
                                 + rhs[index2 + M + 1]
                                 + rhs[index2 + MN]
                                 + rhs[index2 + MN + 1]
                                 + rhs[index2 + MN + M]
                                 + rhs[index2 + MN + M + 1]));
    
    for (i = 1; i < Mhalf - 1; i++)
    {
      index1 = (0 * Nhalf + 0) * Mhalf + i;
      index2 = (2 * 0 * N + 2 * 0) * M + 2 * i;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 + M]
                                   + rhs[index2 + MN]
                                   + rhs[index2 + MN + M])
                            + 0.25 * (rhs[index2 + 1]
                                      + rhs[index2 - 1]
                                      + rhs[index2 + M + 1]
                                      + rhs[index2 + M - 1]
                                      + rhs[index2 + MN + 1]
                                      + rhs[index2 + MN - 1]
                                      + rhs[index2 + MN + M + 1]
                                      + rhs[index2 + MN + M - 1]));
    }
    
    index1 = (0 * Nhalf + 0) * Mhalf + (Mhalf - 1);
    index2 = (2 * 0 * N + 2 * 0) * M + 2 * (Mhalf - 1);
    rhs_coarse[index1] = (0.5 * (rhs[index2]
                                 + rhs[index2 - 1]
                                 + rhs[index2 + M]
                                 + rhs[index2 + M - 1]
                                 + rhs[index2 + MN]
                                 + rhs[index2 + MN - 1]
                                 + rhs[index2 + MN + M]
                                 + rhs[index2 + MN + M - 1]));
    
    for (j = 1; j < Nhalf - 1; j++)
    {
      index1 = (0 * Nhalf + j) * Mhalf + 0;
      index2 = (2 * 0 * N + 2 * j) * M + 2 * 0;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 + 1]
                                   + rhs[index2 + MN]
                                   + rhs[index2 + MN + 1])
                            + 0.25 * (rhs[index2 + M]
                                      + rhs[index2 + M + 1]
                                      + rhs[index2 - M]
                                      + rhs[index2 - M + 1]
                                      + rhs[index2 + MN + M]
                                      + rhs[index2 + MN + M + 1]
                                      + rhs[index2 + MN - M]
                                      + rhs[index2 + MN - M + 1]));
      
      for (i = 1; i < Mhalf - 1; i++)
      {
        index1 = (0 * Nhalf + j) * Mhalf + i;
        index2 = (2 * 0 * N + 2 * j) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + MN])
                              + 0.25 * (rhs[index2 + 1]
                                        + rhs[index2 - 1]
                                        + rhs[index2 + M]
                                        + rhs[index2 - M]
                                        + rhs[index2 + MN + 1]
                                        + rhs[index2 + MN - 1]
                                        + rhs[index2 + MN + M]
                                        + rhs[index2 + MN - M])
                              + 0.125 * (rhs[index2 + M + 1]
                                         + rhs[index2 + M - 1]
                                         + rhs[index2 - M + 1]
                                         + rhs[index2 - M - 1]
                                         + rhs[index2 + MN + M + 1]
                                         + rhs[index2 + MN + M - 1]
                                         + rhs[index2 + MN - M + 1]
                                         + rhs[index2 + MN - M - 1]));
      }
      
      index1 = (0 * Nhalf + j) * Mhalf + (Mhalf - 1);
      index2 = (2 * 0 * N + 2 * j) * M + 2 * (Mhalf - 1);
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 - 1]
                                   + rhs[index2 + MN]
                                   + rhs[index2 + MN - 1])
                            + 0.25 * (rhs[index2 + M]
                                      + rhs[index2 + M - 1]
                                      + rhs[index2 - M]
                                      + rhs[index2 - M - 1]
                                      + rhs[index2 + MN + M]
                                      + rhs[index2 + MN + M - 1]
                                      + rhs[index2 + MN - M]
                                      + rhs[index2 + MN - M - 1]));
    }
    
    index1 = (0 * Nhalf + (Nhalf - 1)) * Mhalf + 0;
    index2 = (2 * 0 * N + 2 * (Nhalf - 1)) * M + 2 * 0;
    rhs_coarse[index1] = (0.5 * (rhs[index2]
                                 + rhs[index2 + 1]
                                 + rhs[index2 - M]
                                 + rhs[index2 - M + 1]
                                 + rhs[index2 + MN]
                                 + rhs[index2 + MN + 1]
                                 + rhs[index2 + MN - M]
                                 + rhs[index2 + MN - M + 1]));
    
    for (i = 1; i < Mhalf - 1; i++)
    {
      index1 = (0 * Nhalf + (Nhalf - 1)) * Mhalf + i;
      index2 = (2 * 0 * N + 2 * (Nhalf - 1)) * M + 2 * i;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 - M]
                                   + rhs[index2 + MN]
                                   + rhs[index2 + MN - M])
                            + 0.25 * (rhs[index2 + 1]
                                      + rhs[index2 - 1]
                                      + rhs[index2 - M + 1]
                                      + rhs[index2 - M - 1]
                                      + rhs[index2 + MN + 1]
                                      + rhs[index2 + MN - 1]
                                      + rhs[index2 + MN - M + 1]
                                      + rhs[index2 + MN - M - 1]));
    }
    
    index1 = (0 * Nhalf + (Nhalf - 1)) * Mhalf + (Mhalf - 1);
    index2 = (2 * 0 * N + 2 * (Nhalf - 1)) * M + 2 * (Mhalf - 1);
    rhs_coarse[index1] = (0.5 * (rhs[index2]
                                 + rhs[index2 - 1]
                                 + rhs[index2 - M]
                                 + rhs[index2 - M - 1]
                                 + rhs[index2 + MN]
                                 + rhs[index2 + MN - 1]
                                 + rhs[index2 + MN - M]
                                 + rhs[index2 + MN - M - 1]));
    
    for (p = 1; p < Phalf - 1; p++)
    {
      index1 = (p * Nhalf + 0) * Mhalf + 0;
      index2 = (2 * p * N + 2 * 0) * M + 2 * 0;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 + 1]
                                   + rhs[index2 + M]
                                   + rhs[index2 + M + 1])
                            + 0.25 * (rhs[index2 + MN]
                                      + rhs[index2 + MN + 1]
                                      + rhs[index2 + MN + M]
                                      + rhs[index2 + MN + M + 1]
                                      + rhs[index2 - MN]
                                      + rhs[index2 - MN + 1]
                                      + rhs[index2 - MN + M]
                                      + rhs[index2 - MN + M + 1]));
      
      for (i = 1; i < Mhalf - 1; i++)
      {
        index1 = (p * Nhalf + 0) * Mhalf + i;
        index2 = (2 * p * N + 2 * 0) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + M])
                              + 0.25 * (rhs[index2 + 1]
                                        + rhs[index2 - 1]
                                        + rhs[index2 + M + 1]
                                        + rhs[index2 + M - 1]
                                        + rhs[index2 + MN]
                                        + rhs[index2 + MN + M]
                                        + rhs[index2 - MN]
                                        + rhs[index2 - MN + M])
                              + 0.125 * (rhs[index2 + MN + 1]
                                         + rhs[index2 + MN - 1]
                                         + rhs[index2 + MN + M + 1]
                                         + rhs[index2 + MN + M - 1]
                                         + rhs[index2 - MN + 1]
                                         + rhs[index2 - MN - 1]
                                         + rhs[index2 - MN + M + 1]
                                         + rhs[index2 - MN + M - 1]));
      }
      
      index1 = (p * Nhalf + 0) * Mhalf + (Mhalf - 1);
      index2 = (2 * p * N + 2 * 0) * M + 2 * (Mhalf - 1);
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 - 1]
                                   + rhs[index2 + M]
                                   + rhs[index2 + M - 1])
                            + 0.25 * (rhs[index2 + MN]
                                      + rhs[index2 + MN - 1]
                                      + rhs[index2 + MN + M]
                                      + rhs[index2 + MN + M - 1]
                                      + rhs[index2 - MN]
                                      + rhs[index2 - MN - 1]
                                      + rhs[index2 - MN + M]
                                      + rhs[index2 - MN + M - 1]));
      
      for (j = 1; j < Nhalf - 1; j++)
      {
        index1 = (p * Nhalf + j) * Mhalf + 0;
        index2 = (2 * p * N + 2 * j) * M + 2 * 0;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 + 1])
                              + 0.25 * (rhs[index2 + M]
                                        + rhs[index2 + M + 1]
                                        + rhs[index2 - M]
                                        + rhs[index2 - M + 1]
                                        + rhs[index2 + MN]
                                        + rhs[index2 + MN + 1]
                                        + rhs[index2 - MN]
                                        + rhs[index2 - MN + 1])
                              + 0.125 * (rhs[index2 + MN + M]
                                         + rhs[index2 + MN + M + 1]
                                         + rhs[index2 + MN - M]
                                         + rhs[index2 + MN - M + 1]
                                         + rhs[index2 - MN + M]
                                         + rhs[index2 - MN + M + 1]
                                         + rhs[index2 - MN - M]
                                         + rhs[index2 - MN - M + 1]));
        
        for (i = 1; i < Mhalf - 1; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          rhs_coarse[index1] = (0.5 * (rhs[index2])
                                + 0.25 * (rhs[index2 + 1]
                                          + rhs[index2 - 1]
                                          + rhs[index2 + M]
                                          + rhs[index2 - M]
                                          + rhs[index2 + MN]
                                          + rhs[index2 - MN])
                                + 0.125 * (rhs[index2 + M + 1]
                                           + rhs[index2 + M - 1]
                                           + rhs[index2 - M + 1]
                                           + rhs[index2 - M - 1]
                                           + rhs[index2 + MN + 1]
                                           + rhs[index2 + MN - 1]
                                           + rhs[index2 + MN + M]
                                           + rhs[index2 + MN - M]
                                           + rhs[index2 - MN + 1]
                                           + rhs[index2 - MN - 1]
                                           + rhs[index2 - MN + M]
                                           + rhs[index2 - MN - M])
                                + 0.0625 * (rhs[index2 + MN + M + 1]
                                            + rhs[index2 + MN + M - 1]
                                            + rhs[index2 + MN - M + 1]
                                            + rhs[index2 + MN - M - 1]
                                            + rhs[index2 - MN + M + 1]
                                            + rhs[index2 - MN + M - 1]
                                            + rhs[index2 - MN - M + 1]
                                            + rhs[index2 - MN - M - 1]));
        }
        
        index1 = (p * Nhalf + j) * Mhalf + (Mhalf - 1);
        index2 = (2 * p * N + 2 * j) * M + 2 * (Mhalf - 1);
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 - 1])
                              + 0.25 * (rhs[index2 + M]
                                        + rhs[index2 + M - 1]
                                        + rhs[index2 - M]
                                        + rhs[index2 - M - 1]
                                        + rhs[index2 + MN]
                                        + rhs[index2 + MN - 1]
                                        + rhs[index2 - MN]
                                        + rhs[index2 - MN - 1])
                              + 0.125 * (rhs[index2 + MN + M]
                                         + rhs[index2 + MN + M - 1]
                                         + rhs[index2 + MN - M]
                                         + rhs[index2 + MN - M - 1]
                                         + rhs[index2 - MN + M]
                                         + rhs[index2 - MN + M - 1]
                                         + rhs[index2 - MN - M]
                                         + rhs[index2 - MN - M - 1]));
      }
      
      index1 = (p * Nhalf + (Nhalf - 1)) * Mhalf + 0;
      index2 = (2 * p * N + 2 * (Nhalf - 1)) * M + 2 * 0;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 + 1]
                                   + rhs[index2 - M]
                                   + rhs[index2 - M + 1])
                            + 0.25 * (rhs[index2 + MN]
                                      + rhs[index2 + MN + 1]
                                      + rhs[index2 + MN - M]
                                      + rhs[index2 + MN - M + 1]
                                      + rhs[index2 - MN]
                                      + rhs[index2 - MN + 1]
                                      + rhs[index2 - MN - M]
                                      + rhs[index2 - MN - M + 1]));
      
      for (i = 1; i < Mhalf - 1; i++)
      {
        index1 = (p * Nhalf + (Nhalf - 1)) * Mhalf + i;
        index2 = (2 * p * N + 2 * (Nhalf - 1)) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 - M])
                              + 0.25 * (rhs[index2 + 1]
                                        + rhs[index2 - 1]
                                        + rhs[index2 - M + 1]
                                        + rhs[index2 - M - 1]
                                        + rhs[index2 + MN]
                                        + rhs[index2 + MN - M]
                                        + rhs[index2 - MN]
                                        + rhs[index2 - MN - M])
                              + 0.125 * (rhs[index2 + MN + 1]
                                         + rhs[index2 + MN - 1]
                                         + rhs[index2 + MN - M + 1]
                                         + rhs[index2 + MN - M - 1]
                                         + rhs[index2 - MN + 1]
                                         + rhs[index2 - MN - 1]
                                         + rhs[index2 - MN - M + 1]
                                         + rhs[index2 - MN - M - 1]));
      }
      
      index1 = (p * Nhalf + (Nhalf - 1)) * Mhalf + (Mhalf - 1);
      index2 = (2 * p * N + 2 * (Nhalf - 1)) * M + 2 * (Mhalf - 1);
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 - 1]
                                   + rhs[index2 - M]
                                   + rhs[index2 - M - 1])
                            + 0.25 * (rhs[index2 + MN]
                                      + rhs[index2 + MN - 1]
                                      + rhs[index2 + MN - M]
                                      + rhs[index2 + MN - M - 1]
                                      + rhs[index2 - MN]
                                      + rhs[index2 - MN - 1]
                                      + rhs[index2 - MN - M]
                                      + rhs[index2 - MN - M - 1]));
    }
    
    index1 = ((Phalf - 1) * Nhalf + 0) * Mhalf + 0;
    index2 = (2 * (Phalf - 1) * N + 2 * 0) * M + 2 * 0;
    rhs_coarse[index1] = (0.5 * (rhs[index2]
                                 + rhs[index2 + 1]
                                 + rhs[index2 + M]
                                 + rhs[index2 + M + 1]
                                 + rhs[index2 - MN]
                                 + rhs[index2 - MN + 1]
                                 + rhs[index2 - MN + M]
                                 + rhs[index2 - MN + M + 1]));
    
    for (i = 1; i < Mhalf - 1; i++)
    {
      index1 = ((Phalf - 1) * Nhalf + 0) * Mhalf + i;
      index2 = (2 * (Phalf - 1) * N + 2 * 0) * M + 2 * i;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 + M]
                                   + rhs[index2 - MN]
                                   + rhs[index2 - MN + M])
                            + 0.25 * (rhs[index2 + 1]
                                      + rhs[index2 - 1]
                                      + rhs[index2 + M + 1]
                                      + rhs[index2 + M - 1]
                                      + rhs[index2 - MN + 1]
                                      + rhs[index2 - MN - 1]
                                      + rhs[index2 - MN + M + 1]
                                      + rhs[index2 - MN + M - 1]));
    }
    
    index1 = ((Phalf - 1) * Nhalf + 0) * Mhalf + (Mhalf - 1);
    index2 = (2 * (Phalf - 1) * N + 2 * 0) * M + 2 * (Mhalf - 1);
    rhs_coarse[index1] = (0.5 * (rhs[index2]
                                 + rhs[index2 - 1]
                                 + rhs[index2 + M]
                                 + rhs[index2 + M - 1]
                                 + rhs[index2 - MN]
                                 + rhs[index2 - MN - 1]
                                 + rhs[index2 - MN + M]
                                 + rhs[index2 - MN + M - 1]));
    
    for (j = 1; j < Nhalf - 1; j++)
    {
      index1 = ((Phalf - 1) * Nhalf + j) * Mhalf + 0;
      index2 = (2 * (Phalf - 1) * N + 2 * j) * M + 2 * 0;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 + 1]
                                   + rhs[index2 - MN]
                                   + rhs[index2 - MN + 1])
                            + 0.25 * (rhs[index2 + M]
                                      + rhs[index2 + M + 1]
                                      + rhs[index2 - M]
                                      + rhs[index2 - M + 1]
                                      + rhs[index2 - MN + M]
                                      + rhs[index2 - MN + M + 1]
                                      + rhs[index2 - MN - M]
                                      + rhs[index2 - MN - M + 1]));
      
      for (i = 1; i < Mhalf - 1; i++)
      {
        index1 = ((Phalf - 1) * Nhalf + j) * Mhalf + i;
        index2 = (2 * (Phalf - 1) * N + 2 * j) * M + 2 * i;
        rhs_coarse[index1] = (0.5 * (rhs[index2]
                                     + rhs[index2 - MN])
                              + 0.25 * (rhs[index2 + 1]
                                        + rhs[index2 - 1]
                                        + rhs[index2 + M]
                                        + rhs[index2 - M]
                                        + rhs[index2 - MN + 1]
                                        + rhs[index2 - MN - 1]
                                        + rhs[index2 - MN + M]
                                        + rhs[index2 - MN - M])
                              + 0.125 * (rhs[index2 + M + 1]
                                         + rhs[index2 + M - 1]
                                         + rhs[index2 - M + 1]
                                         + rhs[index2 - M - 1]
                                         + rhs[index2 - MN + M + 1]
                                         + rhs[index2 - MN + M - 1]
                                         + rhs[index2 - MN - M + 1]
                                         + rhs[index2 - MN - M - 1]));
      }
      
      index1 = ((Phalf - 1) * Nhalf + j) * Mhalf + (Mhalf - 1);
      index2 = (2 * (Phalf - 1) * N + 2 * j) * M + 2 * (Mhalf - 1);
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 - 1]
                                   + rhs[index2 - MN]
                                   + rhs[index2 - MN - 1])
                            + 0.25 * (rhs[index2 + M]
                                      + rhs[index2 + M - 1]
                                      + rhs[index2 - M]
                                      + rhs[index2 - M - 1]
                                      + rhs[index2 - MN + M]
                                      + rhs[index2 - MN + M - 1]
                                      + rhs[index2 - MN - M]
                                      + rhs[index2 - MN - M - 1]));
    }
    
    index1 = ((Phalf - 1) * Nhalf + (Nhalf - 1)) * Mhalf + 0;
    index2 = (2 * (Phalf - 1) * N + 2 * (Nhalf - 1)) * M + 2 * 0;
    rhs_coarse[index1] = (0.5 * (rhs[index2]
                                 + rhs[index2 + 1]
                                 + rhs[index2 - M]
                                 + rhs[index2 - M + 1]
                                 + rhs[index2 - MN]
                                 + rhs[index2 - MN + 1]
                                 + rhs[index2 - MN - M]
                                 + rhs[index2 - MN - M + 1]));
    
    for (i = 1; i < Mhalf - 1; i++)
    {
      index1 = ((Phalf - 1) * Nhalf + (Nhalf - 1)) * Mhalf + i;
      index2 = (2 * (Phalf - 1) * N + 2 * (Nhalf - 1)) * M + 2 * i;
      rhs_coarse[index1] = (0.5 * (rhs[index2]
                                   + rhs[index2 - M]
                                   + rhs[index2 - MN]
                                   + rhs[index2 - MN - M])
                            + 0.25 * (rhs[index2 + 1]
                                      + rhs[index2 - 1]
                                      + rhs[index2 - M + 1]
                                      + rhs[index2 - M - 1]
                                      + rhs[index2 - MN + 1]
                                      + rhs[index2 - MN - 1]
                                      + rhs[index2 - MN - M + 1]
                                      + rhs[index2 - MN - M - 1]));
    }
    
    index1 = ((Phalf - 1) * Nhalf + (Nhalf - 1)) * Mhalf + (Mhalf - 1);
    index2 = (2 * (Phalf - 1) * N + 2 * (Nhalf - 1)) * M + 2 * (Mhalf - 1);
    rhs_coarse[index1] = (0.5 * (rhs[index2]
                                 + rhs[index2 - 1]
                                 + rhs[index2 - M]
                                 + rhs[index2 - M - 1]
                                 + rhs[index2 - MN]
                                 + rhs[index2 - MN - 1]
                                 + rhs[index2 - MN - M]
                                 + rhs[index2 - MN - M - 1]));
  }
}


void
upsample3D(double *rhs, int M, int N, int P,
           double *v, int Mhalf, int Nhalf, int Phalf,
           double *f_out)
{
  int i, j, p;
  int index1;
  int index2;
  int MN = M * N;
  int MNhalf = Mhalf * Nhalf;
  
  if (M % 2 == 0 && N % 2 == 0 && P % 2 == 0)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          if (p == 0)
          {
            if (j == 0)
            {
              if (i == 0)
              {
                f_out[index2] += (v[index1]);
              }
              else
              {
                f_out[index2] += (0.75 * v[index1]
                                  + 0.25 * v[index1 - 1]);
              }
            }
            else
            {
              if (i == 0)
              {
                f_out[index2] += (0.75 * v[index1]
                                  + 0.25 * v[index1 - Mhalf]);
              }
              else
              {
                f_out[index2] += (0.5625 * v[index1]
                                  + 0.1875 * (v[index1 - 1]
                                              + v[index1 - Mhalf])
                                  + 0.0625 * v[index1 - Mhalf - 1]);
              }
            }
          }
          else
          {
            if (j == 0)
            {
              if (i == 0)
              {
                f_out[index2] += (0.75 * v[index1]
                                  + 0.25 * v[index1 - MNhalf]);
              }
              else
              {
                f_out[index2] += (0.5625 * v[index1]
                                  + 0.1875 * (v[index1 - 1]
                                              + v[index1 - MNhalf])
                                  + 0.0625 * v[index1 - MNhalf - 1]);
              }
            }
            else
            {
              if (i == 0)
              {
                f_out[index2] += (0.5625 * v[index1]
                                  + 0.1875 * (v[index1 - Mhalf]
                                              + v[index1 - MNhalf])
                                  + 0.0625 * v[index1 - MNhalf - Mhalf]);
              }
              else
              {
                f_out[index2] += (0.421875 * v[index1]
                                  + 0.140625 * (v[index1 - 1]
                                                + v[index1 - Mhalf]
                                                + v[index1 - MNhalf])
                                  + 0.046875 * (v[index1 - Mhalf - 1]
                                                + v[index1 - MNhalf - 1]
                                                + v[index1 - MNhalf - Mhalf])
                                  + 0.015625 * v[index1 - MNhalf - Mhalf - 1]);
              }
            }
          }
          
          if (p == 0)
          {
            if (j == 0)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + 1] += (0.75 * v[index1]
                                      + 0.25 * v[index1 + 1]);
              }
              else
              {
                f_out[index2 + 1] += (v[index1]);
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + 1] += (0.5625 * v[index1]
                                      + 0.1875 * (v[index1 + 1]
                                                  + v[index1 - Mhalf])
                                      + 0.0625 * v[index1 - Mhalf + 1]);
              }
              else
              {
                f_out[index2 + 1] += (0.75 * v[index1]
                                      + 0.25 * v[index1 - Mhalf]);
              }
            }
          }
          else
          {
            if (j == 0)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + 1] += (0.5625 * v[index1]
                                      + 0.1875 * (v[index1 + 1]
                                                  + v[index1 - MNhalf])
                                      + 0.0625 * v[index1 - MNhalf + 1]);
              }
              else
              {
                f_out[index2 + 1] += (0.75 * v[index1]
                                      + 0.25 * v[index1 - MNhalf]);
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + 1] += (0.421875 * v[index1]
                                      + 0.140625 * (v[index1 + 1]
                                                    + v[index1 - Mhalf]
                                                    + v[index1 - MNhalf])
                                      + 0.046875 * (v[index1 - Mhalf + 1]
                                                    + v[index1 - MNhalf + 1]
                                                    + v[index1 - MNhalf - Mhalf])
                                      + 0.015625 * v[index1 - MNhalf - Mhalf + 1]);
              }
              else
              {
                f_out[index2 + 1] += (0.5625 * v[index1]
                                      + 0.1875 * (v[index1 - Mhalf]
                                                  + v[index1 - MNhalf])
                                      + 0.0625 * v[index1 - MNhalf - Mhalf]);
              }
            }
          }
          
          if (p == 0)
          {
            if (j < Nhalf - 1)
            {
              if (i == 0)
              {
                f_out[index2 + M] += (0.75 * v[index1]
                                      + 0.25 * v[index1 + Mhalf]);
              }
              else
              {
                f_out[index2 + M] += (0.5625 * v[index1]
                                      + 0.1875 * (v[index1 - 1]
                                                  + v[index1 + Mhalf])
                                      + 0.0625 * v[index1 + Mhalf - 1]);
              }
            }
            else
            {
              if (i == 0)
              {
                f_out[index2 + M] += (v[index1]);
              }
              else
              {
                f_out[index2 + M] += (0.75 * v[index1]
                                      + 0.25 * v[index1 - 1]);
              }
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              if (i == 0)
              {
                f_out[index2 + M] += (0.5625 * v[index1]
                                      + 0.1875 * (v[index1 + Mhalf]
                                                  + v[index1 - MNhalf])
                                      + 0.0625 * v[index1 - MNhalf + Mhalf]);
              }
              else
              {
                f_out[index2 + M] += (0.421875 * v[index1]
                                      + 0.140625 * (v[index1 - 1]
                                                    + v[index1 + Mhalf]
                                                    + v[index1 - MNhalf])
                                      + 0.046875 * (v[index1 + Mhalf - 1]
                                                    + v[index1 - MNhalf - 1]
                                                    + v[index1 - MNhalf + Mhalf])
                                      + 0.015625 * v[index1 - MNhalf + Mhalf - 1]);
              }
            }
            else
            {
              if (i == 0)
              {
                f_out[index2 + M] += (0.75 * v[index1]
                                      + 0.25 * v[index1 - MNhalf]);
              }
              else
              {
                f_out[index2 + M] += (0.5625 * v[index1]
                                      + 0.1875 * (v[index1 - 1]
                                                  + v[index1 - MNhalf])
                                      + 0.0625 * v[index1 - MNhalf - 1]);
              }
            }
          }
          
          if (p == 0)
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + M + 1] += (0.5625 * v[index1]
                                          + 0.1875 * (v[index1 + 1]
                                                      + v[index1 + Mhalf])
                                          + 0.0625 * v[index1 + Mhalf + 1]);
              }
              else
              {
                f_out[index2 + M + 1] += (0.75 * v[index1]
                                          + 0.25 * v[index1 + Mhalf]);
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + M + 1] += (0.75 * v[index1]
                                          + 0.25 * v[index1 + 1]);
              }
              else
              {
                f_out[index2 + M + 1] += (v[index1]);
              }
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + M + 1] += (0.421875 * v[index1]
                                          + 0.140625 * (v[index1 + 1]
                                                        + v[index1 + Mhalf]
                                                        + v[index1 - MNhalf])
                                          + 0.046875 * (v[index1 + Mhalf + 1]
                                                        + v[index1 - MNhalf + 1]
                                                        + v[index1 - MNhalf + Mhalf])
                                          + 0.015625 * v[index1 - MNhalf + Mhalf + 1]);
              }
              else
              {
                f_out[index2 + M + 1] += (0.5625 * v[index1]
                                          + 0.1875 * (v[index1 + Mhalf]
                                                      + v[index1 - MNhalf])
                                          + 0.0625 * v[index1 - MNhalf + Mhalf]);
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + M + 1] += (0.5625 * v[index1]
                                          + 0.1875 * (v[index1 + 1]
                                                      + v[index1 - MNhalf])
                                          + 0.0625 * v[index1 - MNhalf + 1]);
              }
              else
              {
                f_out[index2 + M + 1] += (0.75 * v[index1]
                                          + 0.25 * v[index1 - MNhalf]);
              }
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j == 0)
            {
              if (i == 0)
              {
                f_out[index2 + MN] += (0.75 * v[index1]
                                       + 0.25 * v[index1 + MNhalf]);
              }
              else
              {
                f_out[index2 + MN] += (0.5625 * v[index1]
                                       + 0.1875 * (v[index1 - 1]
                                                   + v[index1 + MNhalf])
                                       + 0.0625 * v[index1 + MNhalf - 1]);
              }
            }
            else
            {
              if (i == 0)
              {
                f_out[index2 + MN] += (0.5625 * v[index1]
                                       + 0.1875 * (v[index1 - Mhalf]
                                                   + v[index1 + MNhalf])
                                       + 0.0625 * v[index1 + MNhalf - Mhalf]);
              }
              else
              {
                f_out[index2 + MN] += (0.421875 * v[index1]
                                       + 0.140625 * (v[index1 - 1]
                                                     + v[index1 - Mhalf]
                                                     + v[index1 + MNhalf])
                                       + 0.046875 * (v[index1 - Mhalf - 1]
                                                     + v[index1 + MNhalf - 1]
                                                     + v[index1 + MNhalf - Mhalf])
                                       + 0.015625 * v[index1 + MNhalf - Mhalf - 1]);
              }
            }
          }
          else
          {
            if (j == 0)
            {
              if (i == 0)
              {
                f_out[index2 + MN] += (v[index1]);
              }
              else
              {
                f_out[index2 + MN] += (0.75 * v[index1]
                                       + 0.25 * v[index1 - 1]);
              }
            }
            else
            {
              if (i == 0)
              {
                f_out[index2 + MN] += (0.75 * v[index1]
                                       + 0.25 * v[index1 - Mhalf]);
              }
              else
              {
                f_out[index2 + MN] += (0.5625 * v[index1]
                                       + 0.1875 * (v[index1 - 1]
                                                   + v[index1 - Mhalf])
                                       + 0.0625 * v[index1 - Mhalf - 1]);
              }
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j == 0)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + 1] += (0.5625 * v[index1]
                                           + 0.1875 * (v[index1 + 1]
                                                       + v[index1 + MNhalf])
                                           + 0.0625 * v[index1 + MNhalf + 1]);
              }
              else
              {
                f_out[index2 + MN + 1] += (0.75 * v[index1]
                                           + 0.25 * v[index1 + MNhalf]);
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + 1] += (0.421875 * v[index1]
                                           + 0.140625 * (v[index1 + 1]
                                                         + v[index1 - Mhalf]
                                                         + v[index1 + MNhalf])
                                           + 0.046875 * (v[index1 - Mhalf + 1]
                                                         + v[index1 + MNhalf + 1]
                                                         + v[index1 + MNhalf - Mhalf])
                                           + 0.015625 * v[index1 + MNhalf - Mhalf + 1]);
              }
              else
              {
                f_out[index2 + MN + 1] += (0.5625 * v[index1]
                                           + 0.1875 * (v[index1 - Mhalf]
                                                       + v[index1 + MNhalf])
                                           + 0.0625 * v[index1 + MNhalf - Mhalf]);
              }
            }
          }
          else
          {
            if (j == 0)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + 1] += (0.75 * v[index1]
                                           + 0.25 * v[index1 + 1]);
              }
              else
              {
                f_out[index2 + MN + 1] += (v[index1]);
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + 1] += (0.5625 * v[index1]
                                           + 0.1875 * (v[index1 + 1]
                                                       + v[index1 - Mhalf])
                                           + 0.0625 * v[index1 - Mhalf + 1]);
              }
              else
              {
                f_out[index2 + MN + 1] += (0.75 * v[index1]
                                           + 0.25 * v[index1 - Mhalf]);
              }
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              if (i == 0)
              {
                f_out[index2 + MN + M] += (0.5625 * v[index1]
                                           + 0.1875 * (v[index1 + Mhalf]
                                                       + v[index1 + MNhalf])
                                           + 0.0625 * v[index1 + MNhalf + Mhalf]);
              }
              else
              {
                f_out[index2 + MN + M] += (0.421875 * v[index1]
                                           + 0.140625 * (v[index1 - 1]
                                                         + v[index1 + Mhalf]
                                                         + v[index1 + MNhalf])
                                           + 0.046875 * (v[index1 + Mhalf - 1]
                                                         + v[index1 + MNhalf - 1]
                                                         + v[index1 + MNhalf + Mhalf])
                                           + 0.015625 * v[index1 + MNhalf + Mhalf - 1]);
              }
            }
            else
            {
              if (i == 0)
              {
                f_out[index2 + MN + M] += (0.75 * v[index1]
                                           + 0.25 * v[index1 + MNhalf]);
              }
              else
              {
                f_out[index2 + MN + M] += (0.5625 * v[index1]
                                           + 0.1875 * (v[index1 - 1]
                                                       + v[index1 + MNhalf])
                                           + 0.0625 * v[index1 + MNhalf - 1]);
              }
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              if (i == 0)
              {
                f_out[index2 + MN + M] += (0.75 * v[index1]
                                           + 0.25 * v[index1 + Mhalf]);
              }
              else
              {
                f_out[index2 + MN + M] += (0.5625 * v[index1]
                                           + 0.1875 * (v[index1 - 1]
                                                       + v[index1 + Mhalf])
                                           + 0.0625 * v[index1 + Mhalf - 1]);
              }
            }
            else
            {
              if (i == 0)
              {
                f_out[index2 + MN + M] += (v[index1]);
              }
              else
              {
                f_out[index2 + MN + M] += (0.75 * v[index1]
                                           + 0.25 * v[index1 - 1]);
              }
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.421875 * v[index1]
                                               + 0.140625 * (v[index1 + 1]
                                                             + v[index1 + Mhalf]
                                                             + v[index1 + MNhalf])
                                               + 0.046875 * (v[index1 + Mhalf + 1]
                                                             + v[index1 + MNhalf + 1]
                                                             + v[index1 + MNhalf + Mhalf])
                                               + 0.015625 * v[index1 + MNhalf + Mhalf + 1]);
              }
              else
              {
                f_out[index2 + MN + M + 1] += (0.5625 * v[index1]
                                               + 0.1875 * (v[index1 + Mhalf]
                                                           + v[index1 + MNhalf])
                                               + 0.0625 * v[index1 + MNhalf + Mhalf]);
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.5625 * v[index1]
                                               + 0.1875 * (v[index1 + 1]
                                                           + v[index1 + MNhalf])
                                               + 0.0625 * v[index1 + MNhalf + 1]);
              }
              else
              {
                f_out[index2 + MN + M + 1] += (0.75 * v[index1]
                                               + 0.25 * v[index1 + MNhalf]);
              }
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.5625 * v[index1]
                                               + 0.1875 * (v[index1 + 1]
                                                           + v[index1 + Mhalf])
                                               + 0.0625 * v[index1 + Mhalf + 1]);
              }
              else
              {
                f_out[index2 + MN + M + 1] += (0.75 * v[index1]
                                               + 0.25 * v[index1 + Mhalf]);
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.75 * v[index1]
                                               + 0.25 * v[index1 + 1]);
              }
              else
              {
                f_out[index2 + MN + M + 1] += (v[index1]);
              }
            }
          }
        }
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 0 && P % 2 == 0)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          if (p == 0)
          {
            if (j == 0)
            {
              f_out[index2] += (v[index1]);
            }
            else
            {
              f_out[index2] += (0.75 * v[index1]
                                + 0.25 * v[index1 - Mhalf]);
            }
          }
          else
          {
            if (j == 0)
            {
              f_out[index2] += (0.75 * v[index1]
                                + 0.25 * v[index1 - MNhalf]);
            }
            else
            {
              f_out[index2] += (0.5625 * v[index1]
                                + 0.1875 * (v[index1 - Mhalf]
                                            + v[index1 - MNhalf])
                                + 0.0625 * v[index1 - MNhalf - Mhalf]);
            }
          }
          
          if (p == 0)
          {
            if (j == 0)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + 1] += (0.5 * (v[index1]
                                             + v[index1 + 1]));
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + 1] += (0.375 * (v[index1]
                                               + v[index1 + 1])
                                      + 0.125 * (v[index1 - Mhalf]
                                                 + v[index1 - Mhalf + 1]));
              }
            }
          }
          else
          {
            if (j == 0)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + 1] += (0.375 * (v[index1]
                                               + v[index1 + 1])
                                      + 0.125 * (v[index1 - MNhalf]
                                                 + v[index1 - MNhalf + 1]));
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + 1] += (0.28125 * (v[index1]
                                                 + v[index1 + 1])
                                      + 0.09375 * (v[index1 - Mhalf]
                                                   + v[index1 - Mhalf + 1]
                                                   + v[index1 - MNhalf]
                                                   + v[index1 - MNhalf + 1])
                                      + 0.03125 * (v[index1 - MNhalf - Mhalf]
                                                   + v[index1 - MNhalf - Mhalf + 1]));
              }
            }
          }
          
          if (p == 0)
          {
            if (j < Nhalf - 1)
            {
              f_out[index2 + M] += (0.75 * v[index1]
                                    + 0.25 * v[index1 + Mhalf]);
            }
            else
            {
              f_out[index2 + M] += (v[index1]);
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              f_out[index2 + M] += (0.5625 * v[index1]
                                    + 0.1875 * (v[index1 + Mhalf]
                                                + v[index1 - MNhalf])
                                    + 0.0625 * v[index1 - MNhalf + Mhalf]);
            }
            else
            {
              f_out[index2 + M] += (0.75 * v[index1]
                                    + 0.25 * v[index1 - MNhalf]);
            }
          }
          
          if (p == 0)
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + M + 1] += (0.375 * (v[index1]
                                                   + v[index1 + 1])
                                          + 0.125 * (v[index1 + Mhalf]
                                                     + v[index1 + Mhalf + 1]));
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + M + 1] += (0.5 * (v[index1]
                                                 + v[index1 + 1]));
              }
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + M + 1] += (0.28125 * (v[index1]
                                                     + v[index1 + 1])
                                          + 0.09375 * (v[index1 + Mhalf]
                                                       + v[index1 + Mhalf + 1]
                                                       + v[index1 - MNhalf]
                                                       + v[index1 - MNhalf + 1])
                                          + 0.03125 * (v[index1 - MNhalf + Mhalf]
                                                       + v[index1 - MNhalf + Mhalf + 1]));
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + M + 1] += (0.375 * (v[index1]
                                                   + v[index1 + 1])
                                          + 0.125 * (v[index1 - MNhalf]
                                                     + v[index1 - MNhalf + 1]));
              }
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j == 0)
            {
              f_out[index2 + MN] += (0.75 * v[index1]
                                     + 0.25 * v[index1 + MNhalf]);
            }
            else
            {
              f_out[index2 + MN] += (0.5625 * v[index1]
                                     + 0.1875 * (v[index1 - Mhalf]
                                                 + v[index1 + MNhalf])
                                     + 0.0625 * v[index1 + MNhalf - Mhalf]);
            }
          }
          else
          {
            if (j == 0)
            {
              f_out[index2 + MN] += (v[index1]);
            }
            else
            {
              f_out[index2 + MN] += (0.75 * v[index1]
                                     + 0.25 * v[index1 - Mhalf]);
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j == 0)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + 1] += (0.375 * (v[index1]
                                                    + v[index1 + 1])
                                           + 0.125 * (v[index1 + MNhalf]
                                                      + v[index1 + MNhalf + 1]));
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + 1] += (0.28125 * (v[index1]
                                                      + v[index1 + 1])
                                           + 0.09375 * (v[index1 - Mhalf]
                                                        + v[index1 - Mhalf + 1]
                                                        + v[index1 + MNhalf]
                                                        + v[index1 + MNhalf + 1])
                                           + 0.03125 * (v[index1 + MNhalf - Mhalf]
                                                        + v[index1 + MNhalf - Mhalf + 1]));
              }
            }
          }
          else
          {
            if (j == 0)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + 1] += (0.5 * (v[index1]
                                                  + v[index1 + 1]));
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + 1] += (0.375 * (v[index1]
                                                    + v[index1 + 1])
                                           + 0.125 * (v[index1 - Mhalf]
                                                      + v[index1 - Mhalf + 1]));
              }
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              f_out[index2 + MN + M] += (0.5625 * v[index1]
                                         + 0.1875 * (v[index1 + Mhalf]
                                                     + v[index1 + MNhalf])
                                         + 0.0625 * v[index1 + MNhalf + Mhalf]);
            }
            else
            {
              f_out[index2 + MN + M] += (0.75 * v[index1]
                                         + 0.25 * v[index1 + MNhalf]);
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              f_out[index2 + MN + M] += (0.75 * v[index1]
                                         + 0.25 * v[index1 + Mhalf]);
            }
            else
            {
              f_out[index2 + MN + M] += (v[index1]);
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.28125 * (v[index1]
                                                          + v[index1 + 1])
                                               + 0.09375 * (v[index1 + Mhalf]
                                                            + v[index1 + Mhalf + 1]
                                                            + v[index1 + MNhalf]
                                                            + v[index1 + MNhalf + 1])
                                               + 0.03125 * (v[index1 + MNhalf + Mhalf]
                                                            + v[index1 + MNhalf + Mhalf + 1]));
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.375 * (v[index1]
                                                        + v[index1 + 1])
                                               + 0.125 * (v[index1 + MNhalf]
                                                          + v[index1 + MNhalf + 1]));
              }
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.375 * (v[index1]
                                                        + v[index1 + 1])
                                               + 0.125 * (v[index1 + Mhalf]
                                                          + v[index1 + Mhalf + 1]));
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.5 * (v[index1]
                                                      + v[index1 + 1]));
              }
            }
          }
        }
      }
    }
  }
  
  if (M % 2 == 0 && N % 2 == 1 && P % 2 == 0)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          if (p == 0)
          {
            if (i == 0)
            {
              f_out[index2] += (v[index1]);
            }
            else
            {
              f_out[index2] += (0.75 * v[index1]
                                + 0.25 * v[index1 - 1]);
            }
          }
          else
          {
            if (i == 0)
            {
              f_out[index2] += (0.75 * v[index1]
                                + 0.25 * v[index1 - MNhalf]);
            }
            else
            {
              f_out[index2] += (0.5625 * v[index1]
                                + 0.1875 * (v[index1 - 1]
                                            + v[index1 - MNhalf])
                                + 0.0625 * v[index1 - MNhalf - 1]);
            }
          }
          
          if (p == 0)
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + 1] += (0.75 * v[index1]
                                    + 0.25 * v[index1 + 1]);
            }
            else
            {
              f_out[index2 + 1] += (v[index1]);
            }
          }
          else
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + 1] += (0.5625 * v[index1]
                                    + 0.1875 * (v[index1 + 1]
                                                + v[index1 - MNhalf])
                                    + 0.0625 * v[index1 - MNhalf + 1]);
            }
            else
            {
              f_out[index2 + 1] += (0.75 * v[index1]
                                    + 0.25 * v[index1 - MNhalf]);
            }
          }
          
          if (p == 0)
          {
            if (j < Nhalf - 1)
            {
              if (i == 0)
              {
                f_out[index2 + M] += (0.5 * (v[index1]
                                             + v[index1 + Mhalf]));
              }
              else
              {
                f_out[index2 + M] += (0.375 * (v[index1]
                                               + v[index1 + Mhalf])
                                      + 0.125 * (v[index1 - 1]
                                                 + v[index1 + Mhalf - 1]));
              }
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              if (i == 0)
              {
                f_out[index2 + M] += (0.375 * (v[index1]
                                               + v[index1 + Mhalf])
                                      + 0.125 * (v[index1 - MNhalf]
                                                 + v[index1 - MNhalf + Mhalf]));
              }
              else
              {
                f_out[index2 + M] += (0.28125 * (v[index1]
                                                 + v[index1 + Mhalf])
                                      + 0.09375 * (v[index1 - 1]
                                                   + v[index1 + Mhalf - 1]
                                                   + v[index1 - MNhalf]
                                                   + v[index1 - MNhalf + Mhalf])
                                      + 0.03125 * (v[index1 - MNhalf - 1]
                                                   + v[index1 - MNhalf + Mhalf - 1]));
              }
            }
          }
          
          if (p == 0)
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + M + 1] += (0.375 * (v[index1]
                                                   + v[index1 + Mhalf])
                                          + 0.125 * (v[index1 + 1]
                                                     + v[index1 + Mhalf + 1]));
              }
              else
              {
                f_out[index2 + M + 1] += (0.5 * (v[index1]
                                                 + v[index1 + Mhalf]));
              }
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + M + 1] += (0.28125 * (v[index1]
                                                     + v[index1 + Mhalf])
                                          + 0.09375 * (v[index1 + 1]
                                                       + v[index1 + Mhalf + 1]
                                                       + v[index1 - MNhalf]
                                                       + v[index1 - MNhalf + Mhalf])
                                          + 0.03125 * (v[index1 - MNhalf + 1]
                                                       + v[index1 - MNhalf + Mhalf + 1]));
              }
              else
              {
                f_out[index2 + M + 1] += (0.375 * (v[index1]
                                                   + v[index1 + Mhalf])
                                          + 0.125 * (v[index1 - MNhalf]
                                                     + v[index1 - MNhalf + Mhalf]));
              }
            }
          }
          
          if (p < Phalf - 1)
          {
            if (i == 0)
            {
              f_out[index2 + MN] += (0.75 * v[index1]
                                     + 0.25 * v[index1 + MNhalf]);
            }
            else
            {
              f_out[index2 + MN] += (0.5625 * v[index1]
                                     + 0.1875 * (v[index1 - 1]
                                                 + v[index1 + MNhalf])
                                     + 0.0625 * v[index1 + MNhalf - 1]);
            }
          }
          else
          {
            if (i == 0)
            {
              f_out[index2 + MN] += (v[index1]);
            }
            else
            {
              f_out[index2 + MN] += (0.75 * v[index1]
                                     + 0.25 * v[index1 - 1]);
            }
          }
          
          if (p < Phalf - 1)
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + MN + 1] += (0.5625 * v[index1]
                                         + 0.1875 * (v[index1 + 1]
                                                     + v[index1 + MNhalf])
                                         + 0.0625 * v[index1 + MNhalf + 1]);
            }
            else
            {
              f_out[index2 + MN + 1] += (0.75 * v[index1]
                                         + 0.25 * v[index1 + MNhalf]);
            }
          }
          else
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + MN + 1] += (0.75 * v[index1]
                                         + 0.25 * v[index1 + 1]);
            }
            else
            {
              f_out[index2 + MN + 1] += (v[index1]);
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              if (i == 0)
              {
                f_out[index2 + MN + M] += (0.375 * (v[index1]
                                                    + v[index1 + Mhalf])
                                           + 0.125 * (v[index1 + MNhalf]
                                                      + v[index1 + MNhalf + Mhalf]));
              }
              else
              {
                f_out[index2 + MN + M] += (0.28125 * (v[index1]
                                                      + v[index1 + Mhalf])
                                           + 0.09375 * (v[index1 - 1]
                                                        + v[index1 + Mhalf - 1]
                                                        + v[index1 + MNhalf]
                                                        + v[index1 + MNhalf + Mhalf])
                                           + 0.03125 * (v[index1 + MNhalf - 1]
                                                        + v[index1 + MNhalf + Mhalf - 1]));
              }
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              if (i == 0)
              {
                f_out[index2 + MN + M] += (0.5 * (v[index1]
                                                  + v[index1 + Mhalf]));
              }
              else
              {
                f_out[index2 + MN + M] += (0.375 * (v[index1]
                                                    + v[index1 + Mhalf])
                                           + 0.125 * (v[index1 - 1]
                                                      + v[index1 + Mhalf - 1]));
              }
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.28125 * (v[index1]
                                                          + v[index1 + Mhalf])
                                               + 0.09375 * (v[index1 + 1]
                                                            + v[index1 + Mhalf + 1]
                                                            + v[index1 + MNhalf]
                                                            + v[index1 + MNhalf + Mhalf])
                                               + 0.03125 * (v[index1 + MNhalf + 1]
                                                            + v[index1 + MNhalf + Mhalf + 1]));
              }
              else
              {
                f_out[index2 + MN + M + 1] += (0.375 * (v[index1]
                                                        + v[index1 + Mhalf])
                                               + 0.125 * (v[index1 + MNhalf]
                                                          + v[index1 + MNhalf + Mhalf]));
              }
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.375 * (v[index1]
                                                        + v[index1 + Mhalf])
                                               + 0.125 * (v[index1 + 1]
                                                          + v[index1 + Mhalf + 1]));
              }
              else
              {
                f_out[index2 + MN + M + 1] += (0.5 * (v[index1]
                                                      + v[index1 + Mhalf]));
              }
            }
          }
        }
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 1 && P % 2 == 0)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          if (p == 0)
          {
            f_out[index2] += (v[index1]);
          }
          else
          {
            f_out[index2] += (0.75 * v[index1]
                              + 0.25 * v[index1 - MNhalf]);
          }
          
          if (p == 0)
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + 1] += (0.5 * (v[index1]
                                           + v[index1 + 1]));
            }
          }
          else
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + 1] += (0.375 * (v[index1]
                                             + v[index1 + 1])
                                    + 0.125 * (v[index1 - MNhalf]
                                               + v[index1 - MNhalf + 1]));
            }
          }
          
          if (p == 0)
          {
            if (j < Nhalf - 1)
            {
              f_out[index2 + M] += (0.5 * (v[index1]
                                           + v[index1 + Mhalf]));
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              f_out[index2 + M] += (0.375 * (v[index1]
                                             + v[index1 + Mhalf])
                                    + 0.125 * (v[index1 - MNhalf]
                                               + v[index1 - MNhalf + Mhalf]));
            }
          }
          
          if (p == 0)
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + M + 1] += (0.25 * (v[index1]
                                                  + v[index1 + 1]
                                                  + v[index1 + Mhalf]
                                                  + v[index1 + Mhalf + 1]));
              }
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + M + 1] += (0.1875 * (v[index1]
                                                    + v[index1 + 1]
                                                    + v[index1 + Mhalf]
                                                    + v[index1 + Mhalf + 1])
                                          + 0.0625 * (v[index1 - MNhalf]
                                                      + v[index1 - MNhalf + 1]
                                                      + v[index1 - MNhalf + Mhalf]
                                                      + v[index1 - MNhalf + Mhalf + 1]));
              }
            }
          }
          
          if (p < Phalf - 1)
          {
            f_out[index2 + MN] += (0.75 * v[index1]
                                   + 0.25 * v[index1 + MNhalf]);
          }
          else
          {
            f_out[index2 + MN] += (v[index1]);
          }
          
          if (p < Phalf - 1)
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + MN + 1] += (0.375 * (v[index1]
                                                  + v[index1 + 1])
                                         + 0.125 * (v[index1 + MNhalf]
                                                    + v[index1 + MNhalf + 1]));
            }
          }
          else
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + MN + 1] += (0.5 * (v[index1]
                                                + v[index1 + 1]));
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              f_out[index2 + MN + M] += (0.375 * (v[index1]
                                                  + v[index1 + Mhalf])
                                         + 0.125 * (v[index1 + MNhalf]
                                                    + v[index1 + MNhalf + Mhalf]));
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              f_out[index2 + MN + M] += (0.5 * (v[index1]
                                                + v[index1 + Mhalf]));
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.1875 * (v[index1]
                                                         + v[index1 + 1]
                                                         + v[index1 + Mhalf]
                                                         + v[index1 + Mhalf + 1])
                                               + 0.0625 * (v[index1 + MNhalf]
                                                           + v[index1 + MNhalf + 1]
                                                           + v[index1 + MNhalf + Mhalf]
                                                           + v[index1 + MNhalf + Mhalf + 1]));
              }
            }
          }
          else
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.25 * (v[index1]
                                                       + v[index1 + 1]
                                                       + v[index1 + Mhalf]
                                                       + v[index1 + Mhalf + 1]));
              }
            }
          }
        }
      }
    }
  }
  
  if (M % 2 == 0 && N % 2 == 0 && P % 2 == 1)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          if (j == 0)
          {
            if (i == 0)
            {
              f_out[index2] += (v[index1]);
            }
            else
            {
              f_out[index2] += (0.75 * v[index1]
                                + 0.25 * v[index1 - 1]);
            }
          }
          else
          {
            if (i == 0)
            {
              f_out[index2] += (0.75 * v[index1]
                                + 0.25 * v[index1 - Mhalf]);
            }
            else
            {
              f_out[index2] += (0.5625 * v[index1]
                                + 0.1875 * (v[index1 - 1]
                                            + v[index1 - Mhalf])
                                + 0.0625 * v[index1 - Mhalf - 1]);
            }
          }
          
          if (j == 0)
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + 1] += (0.75 * v[index1]
                                    + 0.25 * v[index1 + 1]);
            }
            else
            {
              f_out[index2 + 1] += (v[index1]);
            }
          }
          else
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + 1] += (0.5625 * v[index1]
                                    + 0.1875 * (v[index1 + 1]
                                                + v[index1 - Mhalf])
                                    + 0.0625 * v[index1 - Mhalf + 1]);
            }
            else
            {
              f_out[index2 + 1] += (0.75 * v[index1]
                                    + 0.25 * v[index1 - Mhalf]);
            }
          }
          
          if (j < Nhalf - 1)
          {
            if (i == 0)
            {
              f_out[index2 + M] += (0.75 * v[index1]
                                    + 0.25 * v[index1 + Mhalf]);
            }
            else
            {
              f_out[index2 + M] += (0.5625 * v[index1]
                                    + 0.1875 * (v[index1 - 1]
                                                + v[index1 + Mhalf])
                                    + 0.0625 * v[index1 + Mhalf - 1]);
            }
          }
          else
          {
            if (i == 0)
            {
              f_out[index2 + M] += (v[index1]);
            }
            else
            {
              f_out[index2 + M] += (0.75 * v[index1]
                                    + 0.25 * v[index1 - 1]);
            }
          }
          
          if (j < Nhalf - 1)
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + M + 1] += (0.5625 * v[index1]
                                        + 0.1875 * (v[index1 + 1]
                                                    + v[index1 + Mhalf])
                                        + 0.0625 * v[index1 + Mhalf + 1]);
            }
            else
            {
              f_out[index2 + M + 1] += (0.75 * v[index1]
                                        + 0.25 * v[index1 + Mhalf]);
            }
          }
          else
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + M + 1] += (0.75 * v[index1]
                                        + 0.25 * v[index1 + 1]);
            }
            else
            {
              f_out[index2 + M + 1] += (v[index1]);
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j == 0)
            {
              if (i == 0)
              {
                f_out[index2 + MN] += (0.5 * (v[index1]
                                              + v[index1 + MNhalf]));
              }
              else
              {
                f_out[index2 + MN] += (0.375 * (v[index1]
                                                + v[index1 + MNhalf])
                                       + 0.125 * (v[index1 - 1]
                                                  + v[index1 + MNhalf - 1]));
              }
            }
            else
            {
              if (i == 0)
              {
                f_out[index2 + MN] += (0.375 * (v[index1]
                                                + v[index1 + MNhalf])
                                       + 0.125 * (v[index1 - Mhalf]
                                                  + v[index1 + MNhalf - Mhalf]));
              }
              else
              {
                f_out[index2 + MN] += (0.28125 * (v[index1]
                                                  + v[index1 + MNhalf])
                                       + 0.09375 * (v[index1 - 1]
                                                    + v[index1 - Mhalf]
                                                    + v[index1 + MNhalf - 1]
                                                    + v[index1 + MNhalf - Mhalf])
                                       + 0.03125 * (v[index1 - Mhalf - 1]
                                                    + v[index1 + MNhalf - Mhalf - 1]));
              }
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j == 0)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + 1] += (0.375 * (v[index1]
                                                    + v[index1 + MNhalf])
                                           + 0.125 * (v[index1 + 1]
                                                      + v[index1 + MNhalf + 1]));
              }
              else
              {
                f_out[index2 + MN + 1] += (0.5 * (v[index1]
                                                  + v[index1 + MNhalf]));
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + 1] += (0.28125 * (v[index1]
                                                      + v[index1 + MNhalf])
                                           + 0.09375 * (v[index1 + 1]
                                                        + v[index1 - Mhalf]
                                                        + v[index1 + MNhalf + 1]
                                                        + v[index1 + MNhalf - Mhalf])
                                           + 0.03125 * (v[index1 - Mhalf + 1]
                                                        + v[index1 + MNhalf - Mhalf + 1]));
              }
              else
              {
                f_out[index2 + MN + 1] += (0.375 * (v[index1]
                                                    + v[index1 + MNhalf])
                                           + 0.125 * (v[index1 - Mhalf]
                                                      + v[index1 + MNhalf - Mhalf]));
              }
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              if (i == 0)
              {
                f_out[index2 + MN + M] += (0.375 * (v[index1]
                                                    + v[index1 + MNhalf])
                                           + 0.125 * (v[index1 + Mhalf]
                                                      + v[index1 + MNhalf + Mhalf]));
              }
              else
              {
                f_out[index2 + MN + M] += (0.28125 * (v[index1]
                                                      + v[index1 + MNhalf])
                                           + 0.09375 * (v[index1 - 1]
                                                        + v[index1 + Mhalf]
                                                        + v[index1 + MNhalf - 1]
                                                        + v[index1 + MNhalf + Mhalf])
                                           + 0.03125 * (v[index1 + Mhalf - 1]
                                                        + v[index1 + MNhalf + Mhalf - 1]));
              }
            }
            else
            {
              if (i == 0)
              {
                f_out[index2 + MN + M] += (0.5 * (v[index1]
                                                  + v[index1 + MNhalf]));
              }
              else
              {
                f_out[index2 + MN + M] += (0.375 * (v[index1]
                                                    + v[index1 + MNhalf])
                                           + 0.125 * (v[index1 - 1]
                                                      + v[index1 + MNhalf - 1]));
              }
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.28125 * (v[index1]
                                                          + v[index1 + MNhalf])
                                               + 0.09375 * (v[index1 + 1]
                                                            + v[index1 + Mhalf]
                                                            + v[index1 + MNhalf + 1]
                                                            + v[index1 + MNhalf + Mhalf])
                                               + 0.03125 * (v[index1 + Mhalf + 1]
                                                            + v[index1 + MNhalf + Mhalf + 1]));
              }
              else
              {
                f_out[index2 + MN + M + 1] += (0.375 * (v[index1]
                                                        + v[index1 + MNhalf])
                                               + 0.125 * (v[index1 + Mhalf]
                                                          + v[index1 + MNhalf + Mhalf]));
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.375 * (v[index1]
                                                        + v[index1 + MNhalf])
                                               + 0.125 * (v[index1 + 1]
                                                          + v[index1 + MNhalf + 1]));
              }
              else
              {
                f_out[index2 + MN + M + 1] += (0.5 * (v[index1]
                                                      + v[index1 + MNhalf]));
              }
            }
          }
        }
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 0 && P % 2 == 1)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          if (j == 0)
          {
            f_out[index2] += (v[index1]);
          }
          else
          {
            f_out[index2] += (0.75 * v[index1]
                              + 0.25 * v[index1 - Mhalf]);
          }
          
          if (j == 0)
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + 1] += (0.5 * (v[index1]
                                           + v[index1 + 1]));
            }
          }
          else
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + 1] += (0.375 * (v[index1]
                                             + v[index1 + 1])
                                    + 0.125 * (v[index1 - Mhalf]
                                               + v[index1 - Mhalf + 1]));
            }
          }
          
          if (j < Nhalf - 1)
          {
            f_out[index2 + M] += (0.75 * v[index1]
                                  + 0.25 * v[index1 + Mhalf]);
          }
          else
          {
            f_out[index2 + M] += (v[index1]);
          }
          
          if (j < Nhalf - 1)
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + M + 1] += (0.375 * (v[index1]
                                                 + v[index1 + 1])
                                        + 0.125 * (v[index1 + Mhalf]
                                                   + v[index1 + Mhalf + 1]));
            }
          }
          else
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + M + 1] += (0.5 * (v[index1]
                                               + v[index1 + 1]));
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j == 0)
            {
              f_out[index2 + MN] += (0.5 * (v[index1]
                                            + v[index1 + MNhalf]));
            }
            else
            {
              f_out[index2 + MN] += (0.375 * (v[index1]
                                              + v[index1 + MNhalf])
                                     + 0.125 * (v[index1 - Mhalf]
                                                + v[index1 + MNhalf - Mhalf]));
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j == 0)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + 1] += (0.25 * (v[index1]
                                                   + v[index1 + 1]
                                                   + v[index1 + MNhalf]
                                                   + v[index1 + MNhalf + 1]));
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + 1] += (0.1875 * (v[index1]
                                                     + v[index1 + 1]
                                                     + v[index1 + MNhalf]
                                                     + v[index1 + MNhalf + 1])
                                           + 0.0625 * (v[index1 - Mhalf]
                                                       + v[index1 - Mhalf + 1]
                                                       + v[index1 + MNhalf - Mhalf]
                                                       + v[index1 + MNhalf - Mhalf + 1]));
              }
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              f_out[index2 + MN + M] += (0.375 * (v[index1]
                                                  + v[index1 + MNhalf])
                                         + 0.125 * (v[index1 + Mhalf]
                                                    + v[index1 + MNhalf + Mhalf]));
            }
            else
            {
              f_out[index2 + MN + M] += (0.5 * (v[index1]
                                                + v[index1 + MNhalf]));
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.1875 * (v[index1]
                                                         + v[index1 + 1]
                                                         + v[index1 + MNhalf]
                                                         + v[index1 + MNhalf + 1])
                                               + 0.0625 * (v[index1 + Mhalf]
                                                           + v[index1 + Mhalf + 1]
                                                           + v[index1 + MNhalf + Mhalf]
                                                           + v[index1 + MNhalf + Mhalf + 1]));
              }
            }
            else
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.25 * (v[index1]
                                                       + v[index1 + 1]
                                                       + v[index1 + MNhalf]
                                                       + v[index1 + MNhalf + 1]));
              }
            }
          }
        }
      }
    }
  }
  
  if (M % 2 == 0 && N % 2 == 1 && P % 2 == 1)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          if (i == 0)
          {
            f_out[index2] += (v[index1]);
          }
          else
          {
            f_out[index2] += (0.75 * v[index1]
                              + 0.25 * v[index1 - 1]);
          }
          
          if (i < Mhalf - 1)
          {
            f_out[index2 + 1] += (0.75 * v[index1]
                                  + 0.25 * v[index1 + 1]);
          }
          else
          {
            f_out[index2 + 1] += (v[index1]);
          }
          
          if (j < Nhalf - 1)
          {
            if (i == 0)
            {
              f_out[index2 + M] += (0.5 * (v[index1]
                                           + v[index1 + Mhalf]));
            }
            else
            {
              f_out[index2 + M] += (0.375 * (v[index1]
                                             + v[index1 + Mhalf])
                                    + 0.125 * (v[index1 - 1]
                                               + v[index1 + Mhalf - 1]));
            }
          }
          
          if (j < Nhalf - 1)
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + M + 1] += (0.375 * (v[index1]
                                                 + v[index1 + Mhalf])
                                        + 0.125 * (v[index1 + 1]
                                                   + v[index1 + Mhalf + 1]));
            }
            else
            {
              f_out[index2 + M + 1] += (0.5 * (v[index1]
                                               + v[index1 + Mhalf]));
            }
          }
          
          if (p < Phalf - 1)
          {
            if (i == 0)
            {
              f_out[index2 + MN] += (0.5 * (v[index1]
                                            + v[index1 + MNhalf]));
            }
            else
            {
              f_out[index2 + MN] += (0.375 * (v[index1]
                                              + v[index1 + MNhalf])
                                     + 0.125 * (v[index1 - 1]
                                                + v[index1 + MNhalf - 1]));
            }
          }
          
          if (p < Phalf - 1)
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + MN + 1] += (0.375 * (v[index1]
                                                  + v[index1 + MNhalf])
                                         + 0.125 * (v[index1 + 1]
                                                    + v[index1 + MNhalf + 1]));
            }
            else
            {
              f_out[index2 + MN + 1] += (0.5 * (v[index1]
                                                + v[index1 + MNhalf]));
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              if (i == 0)
              {
                f_out[index2 + MN + M] += (0.25 * (v[index1]
                                                   + v[index1 + Mhalf]
                                                   + v[index1 + MNhalf]
                                                   + v[index1 + MNhalf + Mhalf]));
              }
              else
              {
                f_out[index2 + MN + M] += (0.1875 * (v[index1]
                                                     + v[index1 + Mhalf]
                                                     + v[index1 + MNhalf]
                                                     + v[index1 + MNhalf + Mhalf])
                                           + 0.0625 * (v[index1 - 1]
                                                       + v[index1 + Mhalf - 1]
                                                       + v[index1 + MNhalf - 1]
                                                       + v[index1 + MNhalf + Mhalf - 1]));
              }
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.1875 * (v[index1]
                                                         + v[index1 + Mhalf]
                                                         + v[index1 + MNhalf]
                                                         + v[index1 + MNhalf + Mhalf])
                                               + 0.0625 * (v[index1 + 1]
                                                           + v[index1 + Mhalf + 1]
                                                           + v[index1 + MNhalf + 1]
                                                           + v[index1 + MNhalf + Mhalf + 1]));
              }
              else
              {
                f_out[index2 + MN + M + 1] += (0.25 * (v[index1]
                                                       + v[index1 + Mhalf]
                                                       + v[index1 + MNhalf]
                                                       + v[index1 + MNhalf + Mhalf]));
              }
            }
          }
        }
      }
    }
  }
  
  if (M % 2 == 1 && N % 2 == 1 && P % 2 == 1)
  {
    for (p = 0; p < Phalf; p++)
    {
      for (j = 0; j < Nhalf; j++)
      {
        for (i = 0; i < Mhalf; i++)
        {
          index1 = (p * Nhalf + j) * Mhalf + i;
          index2 = (2 * p * N + 2 * j) * M + 2 * i;
          f_out[index2] += (v[index1]);
          
          if (i < Mhalf - 1)
          {
            f_out[index2 + 1] += (0.5 * (v[index1]
                                         + v[index1 + 1]));
          }
          
          if (j < Nhalf - 1)
          {
            f_out[index2 + M] += (0.5 * (v[index1]
                                         + v[index1 + Mhalf]));
          }
          
          if (j < Nhalf - 1)
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + M + 1] += (0.25 * (v[index1]
                                                + v[index1 + 1]
                                                + v[index1 + Mhalf]
                                                + v[index1 + Mhalf + 1]));
            }
          }
          
          if (p < Phalf - 1)
          {
            f_out[index2 + MN] += (0.5 * (v[index1]
                                          + v[index1 + MNhalf]));
          }
          
          if (p < Phalf - 1)
          {
            if (i < Mhalf - 1)
            {
              f_out[index2 + MN + 1] += (0.25 * (v[index1]
                                                 + v[index1 + 1]
                                                 + v[index1 + MNhalf]
                                                 + v[index1 + MNhalf + 1]));
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              f_out[index2 + MN + M] += (0.25 * (v[index1]
                                                 + v[index1 + Mhalf]
                                                 + v[index1 + MNhalf]
                                                 + v[index1 + MNhalf + Mhalf]));
            }
          }
          
          if (p < Phalf - 1)
          {
            if (j < Nhalf - 1)
            {
              if (i < Mhalf - 1)
              {
                f_out[index2 + MN + M + 1] += (0.125 * (v[index1]
                                                        + v[index1 + 1]
                                                        + v[index1 + Mhalf]
                                                        + v[index1 + Mhalf + 1]
                                                        + v[index1 + MNhalf]
                                                        + v[index1 + MNhalf + 1]
                                                        + v[index1 + MNhalf + Mhalf]
                                                        + v[index1 + MNhalf + Mhalf + 1]));
              }
            }
          }
        }
      }
    }
  }
}


/* Recursive multigrid function.*/
void
poisson_multigrid3D(double *f, double *d,
		    int n1, int n2, int nm,
		    double *f_out,
		    int M, int N, int P, int *directly_solved)
{
  int i, j, p;
  int k;
  double *r;
  double *r_downsampled;
  double *v;
  int Mhalf;
  int Nhalf;
  int Phalf;
  int MN = M * N;
  
  /* Solve a sufficiently small problem directly. */
  if (M < RECURSION_SIZE_LIMIT
      || N < RECURSION_SIZE_LIMIT
      || P < RECURSION_SIZE_LIMIT)
  {
    solve_directly3D(f, d, f_out, M, N, P);
    *directly_solved = 1;
    return;
  }
  *directly_solved = 0;
  
  /* Initialize solution. */
  memcpy(f_out, f, M * N * P * sizeof(*f_out));
  
  /* Pre-smoothing. */
  for (k = 0; k < n1; k++)
    gauss_seidel3D(f_out, d, M, N, P);
  
  /* Compute residual. */
  r = mxCalloc(M * N * P, sizeof(*r));
  for (p = 0; p < P; p++)
    for (j = 0; j < N; j++)
      for (i = 0; i < M; i++)
      {
	int index = (p * N + j) * M + i;
	double residual = d[index] + 6 * f_out[index];
	if (i == 0)
	  residual -= f_out[index + 1];
	else
	  residual -= f_out[index - 1];
	
	if (i == M - 1)
	  residual -= f_out[index - 1];
	else
	  residual -= f_out[index + 1];
	
	if (j == 0)
	  residual -= f_out[index + M];
	else
	  residual -= f_out[index - M];
	
	if (j == N - 1)
	  residual -= f_out[index - M];
	else
	  residual -= f_out[index + M];
	
	if (p == 0)
	  residual -= f_out[index + MN];
	else
	  residual -= f_out[index - MN];
	
	if (p == P - 1)
	  residual -= f_out[index - MN];
	else
	  residual -= f_out[index + MN];
	
	r[index] = residual;
      }
  
  /* Downsample residual. */
  Mhalf = (M + 1) / 2;
  Nhalf = (N + 1) / 2;
  Phalf = (P + 1) / 2;
  r_downsampled = mxCalloc(Mhalf * Nhalf * Phalf, sizeof(*r_downsampled));
  downsample3D(r, M, N, P, r_downsampled, Mhalf, Nhalf, Phalf);
  
  /* Recurse to compute a correction. */
  v = mxCalloc(Mhalf * Nhalf * Phalf, sizeof(*v));
  for (k = 0; k < nm; k++)
  {
    int directly_solved;
    poisson_multigrid3D(v, r_downsampled, n1, n2, nm, v,
			Mhalf, Nhalf, Phalf, &directly_solved);
    if (directly_solved)
      break;
  }
  
  upsample3D(r, M, N, P, v, Mhalf, Nhalf, Phalf, f_out);
  
  /* Post-smoothing. */
  for (k = 0; k < n2; k++)
    gauss_seidel3D(f_out, d, M, N, P);
  
  mxFree(r);
  mxFree(r_downsampled);
  mxFree(v);
}


/* It is assumed that f_out is initialized to zero when called. */
void
poisson_full_multigrid3D(double *rhs, int number_of_iterations,
			 int M, int N, int P, double *f_out)
{
  double *rhs_downsampled;
  double *f_coarse;
  int k;
  
  /* Unless already coarsest scale, first recurse to coarser scale. */
  if (M >= RECURSION_SIZE_LIMIT
      && N >= RECURSION_SIZE_LIMIT
      && P >= RECURSION_SIZE_LIMIT)
  {
    /* Downsample right hand side. */
    int Mhalf = (M + 1) / 2;
    int Nhalf = (N + 1) / 2;
    int Phalf = (P + 1) / 2;
    rhs_downsampled = mxCalloc(Mhalf * Nhalf * Phalf,
			       sizeof(*rhs_downsampled));
    downsample3D(rhs, M, N, P, rhs_downsampled, Mhalf, Nhalf, Phalf);
    
    f_coarse = mxCalloc(Mhalf * Nhalf * Phalf, sizeof(*f_coarse));
    poisson_full_multigrid3D(rhs_downsampled, number_of_iterations,
			     Mhalf, Nhalf, Phalf, f_coarse);
    /* Upsample the coarse result. */
    upsample3D(rhs, M, N, P, f_coarse, Mhalf, Nhalf, Phalf, f_out);
  }
  
  /* Perform number_of_iterations standard multigrid cycles. */
  for (k = 0; k < number_of_iterations; k++)
  {
    int directly_solved;
    poisson_multigrid3D(f_out, rhs, 2, 2, 2, f_out, M, N, P,
			&directly_solved);
    if (directly_solved)
      break;
  }
}


void
antigradient3D(double *g, double mu, int number_of_iterations,
	       int M, int N, int P, double *f_out)
{
  double *rhs;
  double sum;
  double mean;
  int i, j, p;
  int MN = M * N;
  
  /* Compute right hand side of Poisson problem with Neumann
   * boundary conditions, discretized by finite differences.
   */
  rhs = mxCalloc(M * N * P, sizeof(*rhs));
  for (p = 0; p < P; p++)
    for (j = 0; j < N; j++)
      for (i = 0; i < M; i++)
      {
	int index1 = (p * N + j) * M + i;
	int index2 = index1 + M * N * P;
	int index3 = index1 + 2 * M * N * P;
	double d = 0.0;
	
	if (i == 0)
	  d = g[index1 + 1] + g[index1];
	else if (i == M - 1)
	  d = - g[index1] - g[index1 - 1];
	else
	  d = 0.5 * (g[index1 + 1] - g[index1 - 1]);
	
	if (j == 0)
	  d += g[index2 + M] + g[index2];
	else if (j == N - 1)
	  d += - g[index2] - g[index2 - M];
	else
	  d += 0.5 * (g[index2 + M] - g[index2 - M]);
	
	if (p == 0)
	  d += g[index3 + MN] + g[index3];
	else if (p == P - 1)
	  d += - g[index3] - g[index3 - MN];
	else
	  d += 0.5 * (g[index3 + MN] - g[index3 - MN]);
	
	rhs[index1] = d;
      }
  
  /* Solve the equation system with the full multigrid algorithm.
   * Use W cycles and 2 presmoothing and 2 postsmoothing
   * Gauss-Seidel iterations.
   */
  poisson_full_multigrid3D(rhs, number_of_iterations, M, N, P, f_out);
  
  /* Fix the mean value. */
  sum = 0.0;
  for (i = 0; i < M * N * P; i++)
    sum += f_out[i];
  
  mean = sum / (M * N * P);
  for (i = 0; i < M * N * P; i++)
  {
    f_out[i] -= mean;
    f_out[i] += mu;
  }
}


void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int M, N, P;
  double *g;
  double *f_out;
  double mu;
  int number_of_iterations;
  int dim;
  
  /* Check the input and output arguments. */
  
  /* First we expect an MxNx2 or MxNxPx3 array. */
  dim = mxGetNumberOfDimensions(prhs[0]) - 1;
  if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0])
      || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0])
      || dim < 2 || dim > 3
      || mxGetDimensions(prhs[0])[dim] != dim)
  {
    mexErrMsgTxt("g is expected to be an MxNx2 or MxNxPx3 array.");
  }
  
  M = mxGetDimensions(prhs[0])[0];
  N = mxGetDimensions(prhs[0])[1];
  if (dim > 2)
    P = mxGetDimensions(prhs[0])[2];
  else
    P = 0;
  
  g = mxGetPr(prhs[0]);
  
  /* Next two scalars. */
  if (nrhs < 2)
    mu = 0.0;
  else
  {
    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1])
	|| mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1])
	|| mxGetNumberOfElements(prhs[1]) != 1)
    {
      mexErrMsgTxt("mu is expected to be a scalar.");
    }
    mu = mxGetScalar(prhs[1]);
  }

  if (nrhs < 3)
    number_of_iterations = 2;
  else
  {
    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2])
	|| mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2])
	|| mxGetNumberOfElements(prhs[2]) != 1)
    {
      mexErrMsgTxt("N is expected to be a scalar.");
    }
    number_of_iterations = (int) mxGetScalar(prhs[2]);
    if (number_of_iterations < 0
	|| (double) number_of_iterations != mxGetScalar(prhs[2]))
    {
      mexErrMsgTxt("N is expected to be a positive integer.");
    }
  }
  
  plhs[0] = mxCreateNumericArray(dim, mxGetDimensions(prhs[0]),
				 mxDOUBLE_CLASS, mxREAL);
  f_out = mxGetPr(plhs[0]);
  
  if (dim == 2)
    antigradient2D(g, mu, number_of_iterations, M, N, f_out);
  else
    antigradient3D(g, mu, number_of_iterations, M, N, P, f_out);
}

