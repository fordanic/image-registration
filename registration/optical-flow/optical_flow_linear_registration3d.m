function [transformationMatrix] = ...
    optical_flow_linear_registration3d(moving, fixed, varargin)
% OPTICAL_FLOW_LINEAR_REGISTRATION3D Estimates a transformation matrix using optical flow
%
% INPUT ARGUMENTS
% moving                - Moving iamge
% fixed                 - Fixed image
%
% OUTPUT ARGUMENTS
% transformationMatrix  - Estimated displacement field
%
% OPTIONAL INPUT ARGUMENTS
% 'transformationModel'     - Transformation model for estimating the
%                             displacement field
%                             'translation', 'affine', 'non-rigid' (default)
%
% 'multiModal'              - Set whether to perform multi-modal or
%                             uni-modal image registration
%                             false (default), true
%
% 'numberOfChannels'        - Number of channels to use in when computing
%                             the entropy (based on channel coding). This
%                             is only relevant if multiModal is set to
%                             true.
%                             Default value is 8

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
% translation, affine
transformationModel = 'affine';

% multi-modal
multiModal = false;

% number of channels, only valid for multi-modal registration
numberOfChannels = 8;

% Overwrites default parameter
for k=1:2:length(varargin)
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% Initialize transformation matrix
transformationMatrix = eye(4);

[b_moving(:,:,:,1) b_moving(:,:,:,2) b_moving(:,:,:,3)] = gradient(moving);
[b_fixed(:,:,:,1) b_fixed(:,:,:,2) b_fixed(:,:,:,3)] = gradient(fixed);

if multiModal
    b = b_fixed(:,:,:,2);
    b(:,:,:,2) = b_fixed(:,:,:,1);
    b(:,:,:,3) = b_fixed(:,:,:,3);
    
    [delta_c mask] = estimate_delta_c(fixed,moving,numberOfChannels);
    delta_c(mask ~= 1) = 0;
    mask = repmat(mask,[1 1 3]);
    b(mask ~= 1) = 0;
else
    b = (b_moving(:,:,:,2) + b_fixed(:,:,:,2))/2;
    b(:,:,:,2) = (b_moving(:,:,:,1) + b_fixed(:,:,:,1))/2;
    b(:,:,:,3) = (b_moving(:,:,:,3) + b_fixed(:,:,:,3))/2;
    
    delta_c = fixed - moving;
end

[G, h] = build_G_h_linear3d(b, delta_c, transformationModel);

% Solve the equation system
switch transformationModel
    case 'translation'
        d = G \ h;
        transformationMatrix(1:3,4) = d;
    case {'rigid','affine'}
        p = G \ h;
        transformationMatrix(1:3,1:4) = [1+p(1) p(2) p(3) p(10);...
                                         p(4) 1+p(5) p(6) p(11);...
                                         p(7) p(8) 1+p(9) p(12)];
end