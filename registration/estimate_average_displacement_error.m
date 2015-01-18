function ADE = estimate_average_displacement_error(displacementField1,displacementField2,mask)
% ESTIMATE_AVERAGE_DISPLACEMENT_ERROR Estimates the average displacement error
%
% ADE = estimate_average_displacement_error(displacementField1,displacementField2,mask)
%
% INPUT ARGUMENTS
% displacementField1    - Displacement field 1
% displacementField2    - Displacement field 2
% mask                  - Mask, to remove certain points in the data
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% ADE                   - Average displacement error between the two
%                         provided displacement fields

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

dims = length(displacementField1);
numberOfElements = sum(mask(:) ~= 0);
switch dims
    case 1
        x_error = mask.*(displacementField1{1}-displacementField2{1});
        ADE = sum(sqrt(x_error(:).^2))/numberOfElements;
    case 2
        x_error = mask.*(displacementField1{1}-displacementField2{1});
        y_error = mask.*(displacementField1{2}-displacementField2{2});
        ADE = sum(sqrt(x_error(:).^2 + y_error(:).^2))/numberOfElements;
    case 3
        x_error = mask.*(displacementField1{1}-displacementField2{1});
        y_error = mask.*(displacementField1{2}-displacementField2{2});
        z_error = mask.*(displacementField1{3}-displacementField2{3});
        ADE = sum(sqrt(x_error(:).^2 + y_error(:).^2 + z_error(:).^2))/numberOfElements;
    otherwise
        error('Not implemented for dimensions higher than three.')
end