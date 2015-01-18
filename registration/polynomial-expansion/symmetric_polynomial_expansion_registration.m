function [update] = symmetric_polynomial_expansion_registration(...
    moving, fixed, movingDeformed, fixedDeformed, ...
    multiModal, transformationModel, varargin)
% OPTICAL_FLOW_REGISTRATION Estimates a displacement field using optical flow
%
% [update] = symmetric_optical_flow_registration(...
%     moving, fixed, movingDeformed, fixedDeformed, ...
%     multiModal, transformationModel)
%
% INPUT ARGUMENTS
% moving                        - Moving image
% fixed                         - Fixed image
% movingDeformed                - Deformed moving image
% fixedDeformed                 - Deformed fixed image
% multiModal                    - Multi-modal registration
% transformationModel           - Transformation model for estimating the
%                                 displacement field
%                                 'translation', 'affine', 'non-rigid'
%
% OPTIONAL INPUT ARGUMENTS
% 'signalModel'             - Local signal model to use when computing the
%                             polynomial expansion transformation
%                             'linear' (deafult), 'quadratic'
%
% 'numberOfChannels'        - Number of channels to use in when computing
%                             the entropy (based on channel coding). This
%                             is only relevant if multiModal is set to
%                             true.
%                             Default value is 8
%
% OUTPUT ARGUMENTS
% update
%  displacementUpdateForward    - Estimated update field from fixed to
%                                 moving deformed
%  displacementUpdateBackward   - Estimated update field from moving to
%                                 fixed deformed
% (only if transformation model is set to translation or affine)
%  transformationMatrixForward  - Estimate transformation matrix from fixed 
%                                 to moving deformed
%  transformationMatrixBackward - Estimate transformation matrix from moving 
%                                 to fixed deformed

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

%% Set up default parameters
signalModel = 'linear';
numberOfChannels = 8;

% Overwrites default parameter
for k=1:2:length(varargin)
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

%%
dims = ndims(moving);

if multiModal && strcmp(signalModel,'quadratic')
	error(['Multi-modal image registration using polynomial expansion',...
		' is only available for the linear signal model'])
end

switch dims
    case 2
		[updateForward] = ...
			polynomial_expansion_registration2d(...
			movingDeformed, fixed,...
			'signalModel',signalModel,...
			'transformationModel',transformationModel,...
			'multiModal',multiModal,...
			'numberOfChannels',numberOfChannels);
		[updateBackward] = ...
			polynomial_expansion_registration2d(...
			fixedDeformed, moving,...
			'signalModel',signalModel,...
			'transformationModel',transformationModel,...
			'multiModal',multiModal,...
			'numberOfChannels',numberOfChannels);
    case 3
		[updateForward] = ...
			polynomial_expansion_registration3d(...
			movingDeformed, fixed,...
			'signalModel',signalModel,...
			'transformationModel',transformationModel,...
			'multiModal',multiModal,...
			'numberOfChannels',numberOfChannels);
		[updateBackward] = ...
			polynomial_expansion_registration3d(...
			fixedDeformed, moving,...
			'signalModel',signalModel,...
			'transformationModel',transformationModel,...
			'multiModal',multiModal,...
			'numberOfChannels',numberOfChannels);
end

if strcmp(transformationModel,'non-rigid')
    update.displacementForward = updateForward.displacement;
    update.displacementBackward = updateBackward.displacement;
else
    update.transformationMatrixForward = updateForward.transformationMatrix;
    update.transformationMatrixBackward = updateBackward.transformationMatrix;
end