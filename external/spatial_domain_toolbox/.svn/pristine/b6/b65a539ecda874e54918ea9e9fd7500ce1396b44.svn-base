function f = antigradient(g, mu, n)
    
% ANTIGRADIENT   Reconstruction from gradients in 2D and 3D.
%
% f = antigradient(g) computes the function f which has a gradient
% as close as possible to g. The vector field g must either be an
% MxNx2 array or an MxNxPx3 array.
%
% f = antigradient(g, mu) does the same thing but sets the unrecoverable
% mean of f to mu. When omitted, the mean of f is set to zero.
%
% f = antigradient(g, mu, n) uses n multigrid iterations. The default value
% of 2 iterations is fast but only moderately accurate. The algorithm
% converges quickly, however, so usually 10-15 iterations are sufficient to
% reach full convergence.
%
% This function is based on transforming the inverse gradient problem to a
% PDE - a Poisson equation with Neumann boundary conditions. This equation
% is solved by an implementation of the full multigrid algorithm.
%
% The implementation is most efficient when all sides are of approximately
% the same size. The worst case is when two sides are large and the
% third side much smaller. The reason is that only uniform
% downsampling is implemented and when the smallest side becomes
% smaller than 4, a direct solution is applied. It makes no big
% difference whether the sides are odd or even, although sides which
% are of the form 2^k+1 give somewhat faster convergence.
%
% Author: Gunnar Farnebäck
%         Medical Informatics
%         Linköping University, Sweden
%         gunnar@imt.liu.se

error('ANTIGRADIENT is implemented as a mex-file. It has not been compiled on this platform.')
