function output = deformation(input, displacementField, interpolation)
% DEFORMATION Deforms an image according to the given displacement field
%
% output = deformation(input, displacementField, interpolation)
%
% INPUT ARGUMENTS
% input                 - Data to deform
% displacementField     - Displacement field
% interpolation         - Interpolation method
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% output                - Deformed data

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

if isempty(input)
    output = input;
    return;
end

dims = length(size(input));
sz = size(input);

if (dims == 2)
    [X,Y] = meshgrid(1:sz(2), 1:sz(1));
    output = ba_interp2(input, X+displacementField{1}, ...
        Y+displacementField{2}, interpolation);
elseif (dims == 3)
    [X,Y,Z] = meshgrid(1:sz(2), 1:sz(1), 1:sz(3));
    output = ba_interp3(input, X+displacementField{1}, ...
        Y+displacementField{2}, Z+displacementField{3}, ...
        interpolation);
end
