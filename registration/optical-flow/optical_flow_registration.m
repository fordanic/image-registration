function [update] = optical_flow_registration(...
    moving, fixed, transformationModel, multiModal,varargin)
% OPTICAL_FLOW_REGISTRATION Estimates a displacement field using optical flow
%
% [update] = optical_flow_registration(moving,fixed,transformationModel,multiModal)
%
% INPUT ARGUMENTS
% moving                    - Moving iamge
% fixed                     - Fixed image
% transformationModel       - Transformation model for estimating the
%                             displacement field
%                             'translation', 'affine', 'non-rigid'
% multiModal                - Multi-modal registration
%
% OPTIONAL INPUT ARGUMENTS
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

%% Setup default parameters
% Only valid for multi-modal registration
numberOfChannels = 8;

% Overwrites default parameter
for k=1:2:length(varargin)
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

%%
dims = ndims(moving);

if strcmp(transformationModel,'non-rigid')
    switch dims
        case 2
            [update.displacement, update.certainty] = ...
                demons2d(moving,fixed,...
                'multiModal',multiModal,...
                'numberOfChannels',numberOfChannels);
        case 3
            [update.displacement, update.certainty] = ...
                demons3d(moving,fixed,...
                'multiModal',multiModal,...
                'numberOfChannels',numberOfChannels);
    end
else
    switch dims
        case 2
            [update.transformationMatrix] = ...
                optical_flow_linear_registration2d(moving,fixed,...
                'transformationModel',transformationModel,...
                'multiModal',multiModal,...
                'numberOfChannels',numberOfChannels);
        case 3
            [update.transformationMatrix] = ...
                optical_flow_linear_registration3d(moving,fixed,...
                'transformationModel',transformationModel,...
                'multiModal',multiModal,...
                'numberOfChannels',numberOfChannels);
    end
end
