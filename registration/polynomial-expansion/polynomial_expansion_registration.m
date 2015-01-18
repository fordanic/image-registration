function [update] = polynomial_expansion_registration(...
    moving, fixed, multiModal, transformationModel, varargin)
% POLYNOMIAL_EXPANSION_REGISTRATION Estimates a displacement field using polynomial expansion
%
% INPUT ARGUMENTS
% moving                    - Moving image
% fixed                     - Fixed image
% multiModal                - Set whether to perform multi-modal or
%                             uni-modal image registration
%                             false/true
% transformationModel       - Transformation model for estimating the
%                             displacement field
%                             'translation', 'affine', 'non-rigid'
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
%   displacementUpdate      - Estimated update field
%   certaintyUpdate         - Certainty related to the estimated update field
%   transformationMatrix    - Estimate transformation matrix (only if 

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
		[update] = ...
			polynomial_expansion_registration2d(...
			moving, fixed,...
			'signalModel',signalModel,...
			'transformationModel',transformationModel,...
			'multiModal',multiModal,...
			'numberOfChannels',numberOfChannels);
    case 3
		[update] = ...
			polynomial_expansion_registration3d(...
			moving, fixed,...
			'signalModel',signalModel,...
			'transformationModel',transformationModel,...
			'multiModal',multiModal,...
			'numberOfChannels',numberOfChannels);
end