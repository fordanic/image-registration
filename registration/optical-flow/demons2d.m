function [displacementField, certainty] = demons2d(moving, fixed, varargin)
% DEMONS2D Estimates a displacement field using the original demons scheme
%
% [displacementField, certainty] = demons2d(moving, fixed)
%
% INPUT ARGUMENTS
% moving            - Moving iamge
% fixed             - Fixed image
%
% OPTIONAL INPUT ARGUMENTS
% 'method'                  - Use the gradient of one of the images or both
%                             'fixed', 'moving', 'symmetric'
%
% 'maxDisplacemet'          - Max displacement of the estimated displacement field
%                             2 (default)
%
% 'multiModal'              - Set whether to perform multi-modal or
%                             uni-modal image registration
%                             false (default), true
%                             Setting to true will force method to 'fixed'
%
% 'numberOfChannels'        - Number of channels to use in when computing
%                             the entropy (based on channel coding). This
%                             is only relevant if multiModal is set to
%                             true.
%                             Default value is 8
%
% OUTPUT ARGUMENTS
% displacementField - Estimated displacement field
% certainty         - Certainty map related to the estimated displacement
%                     field


% Copyright (c) 2011
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

% Default parameters
method = 'symmetric';
maxDisplacement = 2;

% multi-modal
multiModal = false;

% Only valid for multi-modal registration
numberOfChannels = 8;

% Overwrites default parameters
for k=1:2:length(varargin),         
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

if (multiModal && ~strcmp(method,'fixed'))
    method = 'fixed';
end

displacementField = cell(2,1);
displacementField{1} = zeros(size(moving));
displacementField{2} = zeros(size(moving));

switch method
    case 'moving'
        % Estimate gradient of moving image
        [movingGradX movingGradY] = gradient(fixed);
        gradX = (movingGradX);
        gradY = (movingGradY);
    case 'fixed'
        % Estimate gradient of fixed image
        [fixedGradX fixedGradY] = gradient(fixed);
        gradX = (fixedGradX);
        gradY = (fixedGradY);
    case 'symmetric'
        % Estimate gradient of fixed image
        [fixedGradX fixedGradY] = gradient(fixed);
        % Estimate gradient of moving image
        [movingGradX movingGradY] = gradient(fixed);
        gradX = 0.5*(fixedGradX + movingGradX);
        gradY = 0.5*(fixedGradY + movingGradY);
end

if multiModal
    numerator = estimate_delta_c(fixed,moving,numberOfChannels);
else
    numerator = (fixed - moving);
end

denominator = numerator.^2 + (gradX.^2 +  gradY.^2) + eps;

displacementField{1} = (numerator.*gradX)./denominator;
displacementField{2} = (numerator.*gradY)./denominator;

currentMaxDisplacement = max(sqrt(displacementField{1}(:).^2 + displacementField{2}(:).^2));

if currentMaxDisplacement > maxDisplacement
    normalizingFactor = currentMaxDisplacement / maxDisplacement;
    displacementField{1} = displacementField{1}/normalizingFactor;
    displacementField{2} = displacementField{2}/normalizingFactor;
end

certainty = sqrt(gradX.^2 +  gradY.^2);