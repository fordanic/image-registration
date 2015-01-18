function [accumulatedField,accumulatedCertainty] = ...
    add_update_field_to_acc_field_with_cert(...
    updateField,accumulatedField,updateCertainty,accumulatedCertainty,...
    accumulationMethod,interpolation)
% ADD_UPDATE_FIELD_TO_ACC_FIELD_WITH_CERT Accumulates the update field to the accumulated field
%
% [accumulatedField] = add_update_field_to_acc_field(...
%   updateField,accumulatedField,,updateCertainty,accumulatedCertainty,...
%   accumulationMethod,interpolation)
%
% INPUT ARGUMENTS
% updateField               - Update field
% accumulatedField          - Accumulated field
% updateCertainty           - Certainty of update field
% accumulatedCertainty      - Certainty of accumulated field
% accumulationMethod        - Accumulation method
% interpolation             - Interpolation to use
%
% OPTIONAL INPUT ARGUMENTS
% N/A
% 
% OUTPUT ARGUMENTS
% accumalatedField 			- Accumulated field
% accumulatedCertainty      - Certainty of accumulated field

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
        certaintySum = updateCertainty + accumulatedCertainty + 1e-16;
        for k = 1 : dims
            accumulatedField{k} = accumulatedField{k} + ...
                updateCertainty.*updateField{k} ./ certaintySum;
        end
        accumulatedCertainty = (accumulatedCertainty.^2 + updateCertainty.^2)./ ...
            certaintySum;
    case 'compositive'
        accumulatedCertainty = deformation(accumulatedCertainty, updateField, interpolation);
        certaintySum = updateCertainty + accumulatedCertainty + 1e-16;
        for k = 1 : dims
            accumulatedField{k} = ...
                deformation(accumulatedField{k}, updateField, interpolation) + ...
                updateCertainty.*updateField{k} ./ certaintySum;
        end
        accumulatedCertainty = (accumulatedCertainty.^2 + updateCertainty.^2)./ ...
            certaintySum;
    case 'diffeomorphic'
        accumulatedCertainty = deformation(accumulatedCertainty, updateField, interpolation);
        certaintySum = updateCertainty + accumulatedCertainty + 1e-16;
        updateField = field_exponentiation(updateField);
        for k = 1 : dims
            accumulatedField{k} = ...
                deformation(accumulatedField{k}, updateField, interpolation) + ...
                updateCertainty.*updateField{k} ./ certaintySum;
        end
        accumulatedCertainty = (accumulatedCertainty.^2 + updateCertainty.^2)./ ...
            certaintySum;
    otherwise
        error('Unknown accumulation method')
end
