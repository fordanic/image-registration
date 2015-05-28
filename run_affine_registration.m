function [deformed, transformationMatrix] = run_affine_registration(...
    moving,fixed,startScale,stopScale,numberOfIterationsPerScale,...
    varargin)
% RUN_AFFINE_REGISTRATION Run affine registration on provided data
%
% run_affine_registration(model,fixed,startScale,stopScale,numberOfIterationsPerScale)
%
% INPUT ARGUMENTS
% moving                        - Moving image
% fixed                         - Fixed image
% startScale                    - Starting scale for registration
% stopScale                     - Stopping scale for registration
% numberOfIterationsPerScale    - Number of iterations per scale
%
% OPTIONAL INPUT ARGUMENTS
% 'method'                      - 'phase'/'optical-flow'/'polynomial-expansion'
%                                 Only relevant for CPU registration
% 'transformationModel'         - 'translation'/'rigid'/'affine'
%                                 Only relevant for CPU registration
% 'symmetric'                   - true/false
% 'rigid'                       - true/false
% 'initialTransformation'       - Initital transformation matrix to apply
% 'numberOfFilters'             - Number of quadrature filters to employ, 3 or 6
%                                 Only relevant for GPU registration
% 'GPU'                         - true/false
% 'useMasks'                    - true/false
% 'movingMask'                  - Certainty mask for moving image
% 'fixedMask'                   - Certainty mask for fixed image
%
% OUTPUT ARGUMENTS
% deformed                      - Transformed moving image
% transformationMatrix          - Transformation matrix

% Copyright (c) 2014 Daniel Forsberg
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

method = 'phase';
transformationModel = 'affine';
symmetric = false;
rigid = false;
initialTransformation = [];
numberOfFilters = 3;
GPU = true;
useMasks = false;
movingMask = [];
fixedMask = [];

% Overwrites default parameters
for k=1:2:length(varargin)
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

%% Registration
if GPU
    %% Set up GPU version
    clear global registrationCUDAAffine;
    
    if ~exist('check_registrationCUDA_affine.m','file') || ...
            ~exist('CUDA_affine_registration3d.m','file')
        error(['Files needed for GPU-based image registration are unavailable. ',...
            'Either they are not on the path or the actual toolboxes for this are missing. ',...
            'The toolboxes can be downloaded at https://bitbucket.org/dforsberg/matlab_cuda ',...
            'and https://bitbucket.org/dforsberg/cuda'])
    end
    
    if ndims(moving) ~= 3 || ndims(fixed) ~= 3
        error('GPU-based image registration is only supported for 3D data.')
    end
    
    global registrationCUDAAffine;
    
    registrationCUDAAffine.moving = moving;
    registrationCUDAAffine.fixed = fixed;
    registrationCUDAAffine.movingMask = movingMask;
    registrationCUDAAffine.fixedMask = fixedMask;
    registrationCUDAAffine.useMasks = useMasks;
    
    % Set registration parameters
    registrationCUDAAffine.dims = ndims(registrationCUDAAffine.moving);
    
    % Number of scales
    registrationCUDAAffine.startScale = startScale;
    registrationCUDAAffine.stopScale = stopScale;
    
    % Iterations per scale
    registrationCUDAAffine.numberOfIterationsPerScale = numberOfIterationsPerScale;
    
    % Symmetric
    registrationCUDAAffine.symmetric = symmetric;
    
    % Rigid
    registrationCUDAAffine.rigid = rigid;
    
    % Initial registration
    registrationCUDAAffine.transformationMatrix = initialTransformation;
    
    
    % Load quadrature filters
    load quadratureFiltersForMorphonRegistration3D
    load quadratureFiltersLinearRegistration3D
    
    if numberOfFilters == 6
        registrationCUDAAffine.quadFilters = f;
    else
        registrationCUDAAffine.quadFilters = {f1, f2, f3};
    end
    
    check_registrationCUDA_affine();
    
    %% Run affine GPU registration
    [deformed, parameters] = CUDA_affine_registration3d(registrationCUDAAffine);
    transformationMatrix = reshape(parameters,[3 4]);
    
    %% Check that transformation matrix is ok, if not reset device and run again
    if isnan(transformationMatrix(1,1)) && numberOfFilters == 6
        warning('Affine CUDA registration failed. Re-running it.')
        registrationCUDAAffine.quadFilters = {f1, f2, f3};
        check_registrationCUDA_affine();
        
        [deformed, parameters] = CUDA_affine_registration3d(registrationCUDAAffine);
        transformationMatrix = reshape(parameters,[3 4]);
        if isnan(transformationMatrix)
            error('Affine CUDA registration failed')
        else
            registrationCUDAAffine.quadFilters = f;
            check_registrationCUDA_affine();
            [deformed, parameters] = CUDA_affine_registration3d(registrationCUDAAffine);
            transformationMatrix = reshape(parameters,[3 4]);
            if isnan(transformationMatrix)
                error('Affine CUDA registration failed')
            end
        end
    end
else
    %% Set up CPU version
    clear global registration;
    
    global registration;
    
    % Set data
    registration.moving = moving;
    registration.fixed = fixed;
    registration.movingMask = movingMask;
    registration.fixedMask = fixedMask;
    registration.useMasks = useMasks;
    registration.dims = ndims(registration.moving);
    
    % Set registration parameters
    registration.method = method;
    registration.transformationModel = transformationModel;
    
    % Number of scales
    registration.startScale = startScale;
    registration.stopScale = stopScale;
    
    % Iterations per scale
    registration.iterationsPerScale = numberOfIterationsPerScale;
    
    % Symmetric
    registration.symmetric = symmetric;
    
    % Rigid
    registration.rigid = rigid;
    
    % Initial registration
    registration.transformationMatrix = initialTransformation;
    
    % Launch registration
    registration_execute();
    
    % Get results
    deformed = registration.deformed;
    transformationMatrix = registration.transformationMatrix;
end
