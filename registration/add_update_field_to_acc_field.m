function [accumulatedField] = add_update_field_to_acc_field(...
    updateField,accumulatedField,accumulationMethod,interpolation)
% ADD_UPDATE_FIELD_TO_ACC_FIELD Accumulates the update field to the accumulated field
%
% [accumulatedField] = add_update_field_to_acc_field(...
%   updateField,accumulatedField,accumulationMethod,interpolation)
%
% INPUT ARGUMENTS
% updateField           - Update field
% accumulatedField      - Accumulated field
% accumulationMethod    - Accumulation method
% interpolation         - Interpolation to use
%
% OPTIONAL INPUT ARGUMENTS
% N/A
% 
% OUTPUT ARGUMENTS
% accumalatedField 		- Accumulated field

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

dims = length(size(updateField{1}));

switch(accumulationMethod)
    case 'sum'
        for k = 1 : dims
            accumulatedField{k} = accumulatedField{k} + updateField{k};                
        end
    case 'compositive'
        for k = 1 : dims
            accumulatedField{k} = deformation(accumulatedField{k}, updateField, interpolation) + updateField{k};                
        end
    case 'diffeomorphic'
        updateField = field_exponentiation(updateField);
        for k = 1 : dims
            accumulatedField{k} = deformation(accumulatedField{k}, updateField, interpolation) + updateField{k};                
        end
    otherwise
        error('Unknown accumulation method')
end
