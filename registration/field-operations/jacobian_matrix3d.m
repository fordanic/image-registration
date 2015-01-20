function [xdx, xdy, xdz, ydx, ydy, ydz, zdx, zdy, zdz] = jacobian_matrix3d(vectorField)
% JACOBIAND_MATRIX3D Estimates the elements of the jacobian matrix
%
% function [xdx xdy xdz ydx ydy ydz zdx zdy zdz] = jacobian_matrix3d(vectorField)
%
% INPUT ARGUMENTS
% vectorField       - Vectorfield to estimate the Jacobian for
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% xdx ...           - Matrix elements of the Jacobian

% Copyright (c) 2012 Daniel Forsberg
% danne.forsberg@outlook.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

[xdx, xdy, xdz] = gradient(vectorField{1});
[ydx, ydy, ydz] = gradient(vectorField{2});
[zdx, zdy, zdz] = gradient(vectorField{2});
