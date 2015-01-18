function [transformationMatrix] = ...
    phase_difference_linear_registration2d(moving, movingCertainty, ...
    fixed, fixedCertainty, varargin)
% PHASE_DIFFERENCE_LINEAR_REGISTRATION2D Estimates a transformation matrix using phase-difference
%
% [transformationMatrix] = ...
%     phase_difference_linear_registration2d(moving, movingCertainty, ...
%     fixed, fixedCertainty)
%
% INPUT ARGUMENTS
% moving                    - Moving image
% movingCertainty           - Certainty mask of moving image
% fixed                     - Fixed image
% fixedCertainty            - Certainty mask of fixed image
%
% OPTIONAL INPUT ARGUMENTS
% 'transformationModel'     - Transformation model for estimating the
%                             displacement field
%                             'translation', 'affine' (default)
%
% OUTPUT ARGUMENTS
% transformationMatrix      - Estimated transformation matrix

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

%% Setup default parameters
% translation, affine, non-rigid
transformationModel = 'affine';

% Overwrites default parameter
for k=1:2:length(varargin)
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

%%
% Initialize transformation matrix
transformationMatrix = eye(3);

% Load quadrature filters
load quadratureFiltersLinearRegistration2D

% Perform quadrature filtering
q21 = imfilter(moving,f1,'same','conv');
q22 = imfilter(moving,f2,'same','conv');

q11 = imfilter(fixed,f1,'same','conv');
q12 = imfilter(fixed,f2,'same','conv');

filterSize = (size(real(f1),1)-1)/2;
q11 = q11(filterSize+1:end-filterSize,filterSize+1:end-filterSize);
q12 = q12(filterSize+1:end-filterSize,filterSize+1:end-filterSize);
q21 = q21(filterSize+1:end-filterSize,filterSize+1:end-filterSize);
q22 = q22(filterSize+1:end-filterSize,filterSize+1:end-filterSize);

% Compute phase-difference, certainties and phase-gradients
dphi(:,:,1) = angle(q11.*conj(q21));
dphi(:,:,2) = angle(q12.*conj(q22));

% Estimate certainties
certainty(:,:,1) = abs(q11.*q21).*((cos(dphi(:,:,1)/2)).^2);
certainty(:,:,2) = abs(q12.*q22).*((cos(dphi(:,:,2)/2)).^2);

% Add certainty masks
if ~isempty(fixedCertainty)
    fixedCertainty = fixedCertainty(...
        filterSize+1:end-filterSize,filterSize+1:end-filterSize);
    certainty = bsxfun(@times,certainty,fixedCertainty);
end
if ~isempty(movingCertainty)
    movingCertainty = movingCertainty(...
        filterSize+1:end-filterSize,filterSize+1:end-filterSize);
    certainty = bsxfun(@times,certainty,movingCertainty);
end

% Estimate gradients of phi
grad_x_dphi_n1 = zeros(size(q11));
grad_x_dphi_n1(:,2:end-1) = ...
    angle(q11(:,3:end).*conj(q11(:,2:end-1)) + ...
    q11(:,2:end-1).*conj(q11(:,1:end-2)) + ...
    q21(:,3:end).*conj(q21(:,2:end-1)) + ...
    q21(:,2:end-1).*conj(q21(:,1:end-2)));

grad_y_dphi_n2 = zeros(size(q11));
grad_y_dphi_n2(2:end-1,:) = ...
    angle(q12(3:end,:).*conj(q12(2:end-1,:)) + ...
    q12(2:end-1,:).*conj(q12(1:end-2,:)) + ...
    q22(3:end,:).*conj(q22(2:end-1,:)) + ...
    q22(2:end-1,:).*conj(q22(1:end-2,:)));

% Save all the phase gradients nicely...
phaseGradient(:,:,1) = grad_x_dphi_n1;
phaseGradient(:,:,2) = grad_y_dphi_n2;

% Build A and h for A p = h
[A, h] = build_A_h_2d(dphi, certainty, phaseGradient, transformationModel);

% Solve the equation system
switch transformationModel
    case 'translation'
        d = A \ h;
        transformationMatrix(1:2,3) = d;
    case {'rigid','affine'}
        p = A \ h;
        transformationMatrix(1:2,1:3) = [1+p(1) p(2) p(5);...
            p(3) 1+p(4) p(6)];
end
