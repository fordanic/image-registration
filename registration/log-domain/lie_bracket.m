function [vectorField] = lie_bracket(vectorField,updateVectorField)
% LIE_BRACKET Computes the Lie bracket of vector fields
%
% [vectorField] = lie_bracket(vectorField,updateVectorField)
%
% INPUT ARGUMENTS
% vectorField           - First vector field
% updateVectorField     - Second vector field, assumed to be small
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% vectorField 			- Final vector field

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

dims = ndims(vectorField{1});
if dims == 2
    [v_xdx, v_xdy, v_ydx, v_ydy] = jacobian_matrix2d(vectorField);
    [u_xdx, u_xdy, u_ydx, u_ydy] = jacobian_matrix2d(updateVectorField);
    temp{1} = (v_xdx.*updateVectorField{1} + v_xdy.*updateVectorField{2}) - ...
        (u_xdx.*vectorField{1} + u_xdy.*vectorField{2});
    temp{2} = (v_ydx.*updateVectorField{1} + v_ydy.*updateVectorField{2}) - ...
        (u_ydx.*vectorField{1} + u_ydy.*vectorField{2});
    vectorField = temp;
else
    [v_xdx, v_xdy, v_xdz, v_ydx, v_ydy, v_ydz, v_zdx, v_zdy, v_zdz] = jacobian_matrix3d(vectorField);
    [u_xdx, u_xdy, u_xdz, u_ydx, u_ydy, u_ydz, u_zdx, u_zdy, u_zdz] = jacobian_matrix3d(updateVectorField);
    temp{1} = (v_xdx.*updateVectorField{1} + v_xdy.*updateVectorField{2} + v_xdz.*updateVectorField{3}) - ...
        (u_xdx.*vectorField{1} + u_xdy.*vectorField{2} + u_xdz.*vectorField{3});
    temp{2} = (v_ydx.*updateVectorField{1} + v_ydy.*updateVectorField{2} + v_ydz.*updateVectorField{3}) - ...
        (u_ydx.*vectorField{1} + u_ydy.*vectorField{2} + u_ydz.*vectorField{3});
    temp{3} = (v_zdx.*updateVectorField{1} + v_zdy.*updateVectorField{2} + v_zdz.*updateVectorField{3}) - ...
        (u_zdx.*vectorField{1} + u_zdy.*vectorField{2} + u_zdz.*vectorField{3});
    vectorField = temp;
end