function [deformed, displacementField, displacementFieldInverse] = ...
    run_non_rigid_registration(...
    moving, fixed, method, accumulationMethod, ...
    startScale, stopScale, numberOfIterationsPerScale, ...
    fluidRegularization, fluidRegularizationData, ...
    elasticRegularization, elasticRegularizationData, ...
    varargin)
% RUN_NON_RIGID_REGISTRATION Run non-rigid registration on provided data
%
% run_non_rigid_registration(moving, fixed, method, accumulationMethod, ...
% 	startScale, stopScale, numberOfIterationsPerScale, ...
% 	fluidRegularization, fluidRegularizationData, ...
%   elasticRegularization, elasticRegularizationData)
%
% INPUT ARGUMENTS
% moving                        - Moving image
% fixed                         - Fixed image
% method                        - 'morphon' or 'demons'
% accumulationMethod            - 'add'/'compositive'/'diffeomorphic'
% startScale                    - Start scale for registration
% stopScale                     - Stop scale for registration
% numberOfIterationsPerScale    - Number of iterations per scale
% fluidRegularization           - Regularization of update displacement field
% fluidRegularizationData       - Sigma for regularization
% elasticRegularization         - Regularization of accumulated displacement field
% elasticRegularizationData     - Sigma for regularization
%
% OPTIONAL INPUT ARGUMENTS
% 'symmetric'                   - true/false
% 'inverse'                     - true/false
% 'GPU'                         - true/false
%
% OUTPUT ARGUMENTS
% deformed                      - Deformed moving image
% displacementField             - Estimated displacement field
% displacementFieldInverse      - Inverse displacement field

% Copyright (c) 2015 Daniel Forsberg
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

%% Default parameters

symmetric = false;
inverse = false;
GPU = true;

% Overwrites default parameters
for k=1:2:length(varargin)
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

%% Registration
if GPU
    %% Set GPU version
    clear global registrationCUDA
    
    if ~exist('check_registrationCUDA_affine.m','file') || ...
            ~exist('CUDA_non_rigid_registration3d.m','file')
        error(['Files needed for GPU-based image registration are unavailable. ',...
            'Either they are not on the path or the actual toolboxes for this are missing. ',...
            'The toolboxes can be downloaded at https://bitbucket.org/dforsberg/matlab_cuda ',...
            'and https://bitbucket.org/dforsberg/cuda'])
    end
    
    if ndims(moving) ~= 3 || ndims(fixed) ~= 3
        error('GPU-based image registration is only supported for 3D data.')
    end
    
    global registrationCUDA;
    
    registrationCUDA.moving = moving;
    registrationCUDA.fixed = fixed;
    
    registrationCUDA.dims = ndims(registrationCUDA.fixed);
    
    registrationCUDA.method = method;
    
    % add, compositive, diffeomorphic
    registrationCUDA.accumulation = accumulationMethod;
    
    % Number of scales
    registrationCUDA.startScale = startScale;
    registrationCUDA.stopScale = stopScale;
    
    % Iterations per scale
    registrationCUDA.numberOfIterationsPerScale = numberOfIterationsPerScale;
    
    % Symmetric
    registrationCUDA.symmetric = symmetric;
    
    % Inverse
    registrationCUDA.inverse = inverse;
    
    % Regularization parameters
    registrationCUDA.fluidRegularization = fluidRegularization;
    registrationCUDA.fluidRegularizationData = fluidRegularizationData;
    registrationCUDA.elasticRegularization = elasticRegularization;
    registrationCUDA.elasticRegularizationData = elasticRegularizationData;
    
    % Morphon specific parameters
    load quadratureFiltersForMorphonRegistration3D
    filterDirection = cell2mat(filterDirection);
    registrationCUDA.filterDirectionX = filterDirection(2:3:end);
    registrationCUDA.filterDirectionY = filterDirection(1:3:end);
    registrationCUDA.filterDirectionZ = filterDirection(3:3:end);
    registrationCUDA.quadFilters = f;
    
    check_registrationCUDA();
    
    %% Run non-rigid GPU registration
    displacementField = cell(3,1);
    displacementFieldInverse = cell(3,1);
    
    if inverse
        [deformed, displacementField{1}, displacementField{2}, displacementField{3}, ...
            displacementFieldInverse{1}, displacementFieldInverse{2}, displacementFieldInverse{3}] = ...
            CUDA_non_rigid_registration3d(registrationCUDA);
    else
        [deformed, displacementField{1}, displacementField{2}, displacementField{3}] = ...
            CUDA_non_rigid_registration3d(registrationCUDA);
    end
else
    %% Set CPU version
    clear global registration
    
    global registration;
    
    % Set data
    registration.moving = moving;
    registration.fixed = fixed;
    registration.dims = ndims(registration.fixed);
    
    % Set registration parameters
    if strcmp(method,'morphon')
        registration.method = 'phase';
    elseif strcmp(method,'demons')
        registration.method = 'optical-flow';
    end
    registration.transformationModel = 'non-rigid';
    
    % add, compositive, diffeomorphic
    registration.accumulationMethod = accumulationMethod;
    
    % Number of scales
    registration.startScale = startScale;
    registration.stopScale = stopScale;
    
    % Iterations per scale
    registration.iterationsPerScale = numberOfIterationsPerScale;
    
    % Symmetric
    registration.symmetric = symmetric;
    
    % Inverse
    registration.inverse = inverse;
    
    % Regularization parameters
    registration.fluidRegularization = fluidRegularization;
    registration.fluidRegularizationData = fluidRegularizationData;
    registration.elasticRegularization = elasticRegularization;
    registration.elasticRegularizationData = elasticRegularizationData;
    
    % Launch registration
    registration_execute();
    
    % Get results
    deformed = registration.deformed;
    displacementField = registration.displacementField;
    if inverse
        displacementFieldInverse = registration.inverseDisplacementField;
    else
        displacementFieldInverse = cell(3,1);
    end
end