function [displacementField] = field_regularization(displacementField, certainty,...
    applyCertainty,sigma)
% FIELD_REGULARIZATION Regularizes the displacement field according to certainty and sigma
%
% [displacementField] = field_regularization(displacementField, certainty, applyCertainty,sigma)
%
% INPUT ARGUMENTS
% displacementField     - Displacement field to regularize
% certainty             - Certainty related to the displacement field
% applyCertainty        - Apply the certainty in the regularization process
% sigma                 - Sigma of the Gaussian kernel applied in the
%                         regularization
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% displacementField     - Regularized displacement field

% Copyright (C) 2012  Daniel Forsberg
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

dims = ndims(displacementField{1});

if applyCertainty
    for k = 1 : dims
        displacementField{k} = normgauss_smoothing(...
            displacementField{k},certainty,sigma);
    end
else
    for k = 1 : dims
        displacementField{k} = gauss_smoothing(...
            displacementField{k},sigma);
    end
end