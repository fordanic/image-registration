function [vectorField] = backward(vectorField)
% BACKWARD Flips the direction of a given vector field
%
% [vectorField] = backward(vectorField)
%
% INPUT ARGUMENTS
% vectorField 	- Vector field to flip
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% vectorField 	- Flipped vector field

% Copyright (c) 2012 Daniel Forsberg
% danne.forsberg@outlook.com

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

for k = 1 : dims
    vectorField{k} = - vectorField{k};
end