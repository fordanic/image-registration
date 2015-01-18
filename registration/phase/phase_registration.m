function [update] = phase_registration(...
    moving, movingCertainty, fixed, fixedCertainty, ...
    transformationModel, currentScale)
% PHASE_REGISTRATION Estimates a displacement field using phase-difference
%
% [update] = phase_registration(...
%     moving, movingCertainty, fixed, fixedCertainty, ...
%     transformationModel, currentScale)
%
% INPUT ARGUMENTS
% moving                        - Moving iamge
% fixed                         - Fixed image
% movingDeformed                - Deformed moving iamge
% fixedDeformed                 - Deformed fixed image
% movingCertainty               - Certainty mask moving iamge
% fixedCertainty                - Certainty mask fixed image
% movingDeformedCertainty       - Certainty mask deformed moving iamge
% fixedDeformedCertainty        - Certainty mask deformed fixed image
% transformationModel           - Transformation model for estimating the
%                                 displacement field
%                                 'translation', 'affine', 'non-rigid'
% currentScale                  - Current scale
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% update
%   displacementUpdate      - Estimated update field
%   certaintyUpdate         - Certainty related to the estimated update field
%   transformationMatrix    - Estimate transformation matrix (only if 
%                             transformation model is set to translation or affine)

% Copyright (c) 2012
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

dims = ndims(moving);

if strcmp(transformationModel,'non-rigid')
    switch dims
        case 2
            [update.displacement, update.certainty] = ...
                morphon2d(moving, movingCertainty, ...
                fixed, fixedCertainty, ...
                currentScale);
        case 3
            [update.displacement, update.certainty] = ...
                morphon3d(moving, movingCertainty, ...
                fixed, fixedCertainty, ...
                currentScale);
    end
else
    switch dims
        case 2
            [update.transformationMatrix] = ...
                phase_difference_linear_registration2d(...
                moving, movingCertainty, ...
                fixed, fixedCertainty, ...
                'transformationModel',transformationModel);
        case 3
            [update.transformationMatrix] = ...
                phase_difference_linear_registration3d(...
                moving, movingCertainty, ...
                fixed, fixedCertainty, ...
                'transformationModel',transformationModel);
    end
end