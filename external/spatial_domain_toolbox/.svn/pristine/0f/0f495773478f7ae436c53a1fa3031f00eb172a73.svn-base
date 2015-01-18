function t = tensorprodc(varargin)

% Do the heavy lifting for tensorprod.m. You can call this function directly
% but you had better know exactly what you do since it's implemented as a C
% mex function with no error checking whatsoever. If you make a mistake
% the whole Matlab process might crash.
% 
% If you do call it, this is the expected format:
%
% t = tensorprodc(sizes, output_indices, arrays..., indices...)
%
% where 'sizes' is a vector (column or row) of the size along each index,
% 'output_indices' is a vector (column or row) with the numbers of the
% indices used in the output, 'arrays...' are the N arrays to be multiplied,
% and 'indices...' are vectors (column only) of the corresponding index
% numbers.
%
% For example,
%
% tensorprod(A, 'ijk', B, 'jmk', C, 'jl')
%
% where A has size (I,J,K), B is (J,M,K), and C is (J,L), is translated into
% the call
%
% tensorprodc([I J K M L], [1;4;5], A, B, C, [1;2;3], [2;4;3], [2;5])
%
% Author: Gunnar Farnebäck
%         Medical Informatics
%         Linköping University, Sweden
%         gunnar@imt.liu.se

error('TENSORPRODC is implemented as a mex-file. It has not been compiled on this platform.')
