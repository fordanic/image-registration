function varargout = neighborhoodloop(varargin)
% NEIGHBORHOODLOOP Implicit loop over neighborhoods.
%
% Example:
% Let f be a matrix. Then
% h = neighborhoodloop(2, [3 5], f, ones(size(f)), @(f,c) median(f(c>0)));
% implements a median filter with 3x5 kernel and proper median computations
% near the borders.
%
% This can be generalized to
% [y1, y2, ...] = neighborhoodloop(N, size, x1, x2, ..., func, a1, a2, ...);
%
% where y1, y2, ... is an arbitrary number of output variables,
% the loops are over the first N dimensions in the arbitrary number
% of arrays x1, x2, ..., and a1, a2, ... is an arbitrary number of
% additional parameters. The called function is given by func,
% which is the name of either a builtin function, an m-file function,
% or a mex function. Alternatively func can be a function handle or
% an inline function. It is supposed that the first argument after
% size which is either a string or a scalar is func. The size
% parameter determines how large neighborhoods are passed to func. It
% can be either a scalar or a vector of length N, all numbers must be
% odd integers. Neighborhood values outside the border are set to
% zero.
%
% See also ARRAYLOOP.
%
% Author: Gunnar Farnebäck
%         Medical Informatics
%         Linköping University, Sweden
%         gunnar@imt.liu.se

error('NEIGHBORHOODLOOP is implemented as a mex-file. It has not been compiled on this platform.')
