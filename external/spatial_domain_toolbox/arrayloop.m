function varargout = arrayloop(varargin)
% ARRAYLOOP Implicit loop over some of the array dimensions.
%
% Example:
% [x, dx] = arrayloop(2, a, b, 'lscov', v);
%
% is equivalent to, but differently implemented than,
%
% for i = 1:size(a, 1)
%     for j = 1:size(a, 2)
% 	   [x(i, j, :) dx(i, j, :)] = ...
% 		lscov(squeeze(a(i, j, :, :)), ...
% 		      squeeze(b(i, j, :)), ...
% 		      v);
%     end
% end
%
% This can be generalized to
% [y1, y2, ...] = arrayloop(N, x1, x2, ..., func, a1, a2, ...);
%
% where y1, y2, ... is an arbitrary number of output variables,
% the loops are over the first N dimensions in the arbitrary number
% of arrays x1, x2, ..., and a1, a2, ... is an arbitrary number of
% additional parameters. The called function is given by func,
% which is the name of either a builtin function, an m-file function,
% or a mex function. Alternatively func can be a function handle or
% an inline function. It is supposed that the first argument after N
% which is either a string or a scalar is func.
%
% See also NEIGHBORHOODLOOP.
%
% Author: Gunnar Farnebäck
%         Medical Informatics
%         Linköping University, Sweden
%         gunnar@imt.liu.se

error('ARRAYLOOP is implemented as a mex-file. It has not been compiled on this platform.')
