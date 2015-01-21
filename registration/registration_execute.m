function registration_execute()
% REGISTRATION_EXECUTE Performs image registration on the data in the struct 'registration'
%
% registration_execute()
%
% Note that the struct 'registration' is a global variable
%
% struct registration
% fixed                             - Fixed image
% moving                            - Moving image
% method                            - Basic concept of the applied method
%                                     'polynomial-expansion', 'optical-flow', 'phase'
% transformationModel               - Transformation model
%                                     'translation', 'rigid', 'affine', (parametric)
%                                     'non-rigid' (non-parametric)
% multiLevelMode                    - Type of multi level mode
%                                     'scale' or 'resolution'
% multiModal                        - Multi-modal image registration, true/false
%                                     Note that this feature is very unstable
%                                     and has not been fully evaluated. Please
%                                     try it but if you do want to perform
%                                     multi-modal image registration, please use
%                                     another toolbox.
% symmetric                         - Symmetric registration, true/false
%                                     If symmetric and non-rigid then
%                                     log-domain will be set to true
% startScale                        - Start scale
% stopScale                         - Stop scale <= start scale
% iterationsPerScale                - Number of iterations per scale
%                                     [# iteration stop scale, ... # iterations start scale]
%
% Only relevant for transformationModel set to parametric registration
% initialTransformationMatrix       - Initial transformation matrix to start
%                                     from
%
% Only relevant for multi-modal and optical flow or polynomial expansion
% numberOfChannels                  - Number of channels to use in channel
%                                     coding when computing the distributions
%
% Only relevant for polynomial-expansion
% signalModel                       - Polynomial signal model used in
%                                     polynomial expansion
%                                     'linear', 'quadratic'
%
% Only relevent for phase
% fixedCertainty                    - Certainty mask associated with the
%                                     fixed image
% movingCertainty                   - Certainty mask associated with the
%                                     moving image
%
% Only relevant for non-rigid registration
% fluidRegularization               - Regularization of the update field, true/false
% fluidRegularizationData           - Sigma to apply per scale
%                                     [# sigma stop scale, ... # sigma start scale]
% elasticRegularization             - Regularization of the accumulated field, true/false
% elasticRegularizationData         - Sigma to apply per scale
%                                     [# sigma stop scale, ... # sigma start scale]
% accumulationMethod                - Method for accumulating the update
%                                     field to the accumulated field
%                                     'sum', 'compositive', 'diffeomorphic'
% applyCertainty                    - Apply certainty when regularizing the
%                                     displacement fields and accumulating
%                                     the displacement fields
%                                     Not implemented yet, i.e. do not use.
% logDomain                         - Perform the registration in the
%                                     log-domain, true/false
% BCHmode                           - Baker-Campbell-Hausdorff
%                                     approximation to use for velocity
%                                     field accumulation in the log-domain
%
% evaluate                          - Evalute various registration metrics during
%                                     the registration process, true/false
% display                           - Display intermediate resulta, true/false
% displayFinal                      - Display final results, true/false
% gamma                             - Applied gamma for visualization, [0,1]
%
% For more details on the parameters and basically describing the whole
% framework used here, please see my thesis "Robust Image Registration for 
% Improved Clinical Efficiency: Using Local Structure Analysis and Model-Based 
% Processing" available at:
% http://urn.kb.se/resolve?urn=urn:nbn:se:liu:diva-91116.

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

%% Set up data
% Get global parameter
global registration;

% Check parameter values
check_registration();

if registration.recordMovie
    vidObj = VideoWriter(registration.movieName);
    vidObj.FrameRate = 2;
    open(vidObj);
end

%% Transform to tensor magnitude
if registration.useTensorMagnitudeRepresentation
    if registration.displayFinal
        fprintf('Transform data to tensor magnitude representation\n');
    end
    switch registration.dims
        case 2
            [t11, t12, t22] = compute_structure_tensor2d(...
                registration.fixed,'average',false,'normalize',true);
            registration.fixed = sqrt(t11.^2 + 2*t12.^2 + t22.^2);
            [t11, t12, t22] = compute_structure_tensor2d(...
                registration.moving,'average',false,'normalize',true);
            registration.moving = sqrt(t11.^2 + 2*t12.^2 + t22.^2);
        case 3
            [t11, t12, t13, t22, t23, t33] = compute_structure_tensor3d(...
                registration.fixed,'average',false,'normalize',true);
            registration.fixed = sqrt(t11.^2 + 2*t12.^2 + 2*t13.^2 + ...
                t22.^2 + 2*t23.^2 + t33.^2);
            [t11, t12, t13, t22, t23, t33] = compute_structure_tensor3d(...
                registration.moving,'average',false,'normalize',true);
            registration.moving = sqrt(t11.^2 + 2*t12.^2 + 2*t13.^2 + ...
                t22.^2 + 2*t23.^2 + t33.^2);
    end
end

%% Create multi-level/scale representation
registration.originalSize = size(registration.fixed);
switch registration.multiLevelMode
    case 'scale'
        if registration.displayFinal
            fprintf('Create multi-scale representation of the data\n');
        end
        registration.fixed = create_multi_scale_representation(...
            registration.fixed, registration.startScale + 1);
        registration.moving = create_multi_scale_representation(...
            registration.moving, registration.startScale + 1);
        if strcmp(registration.method,'phase') && ...
                ~isempty(registration.fixedCertainty)
            registration.fixedCertainty = create_multi_scale_representation(...
                registration.fixedCertainty, registration.startScale + 1);
        else
            registration.fixedCertainty = cell(1,registration.startScale + 1);
        end
        if strcmp(registration.method,'phase') && ...
                ~isempty(registration.movingCertainty)
            registration.movingCertainty = create_multi_scale_representation(...
                registration.movingCertainty, registration.startScale + 1);
        else
            registration.movingCertainty = cell(1,registration.startScale + 1);
        end
    case 'resolution'
        if registration.displayFinal
            fprintf('Create multi-resolution representation of the data\n');
        end
        registration.fixed = create_multi_level_representation(...
            registration.fixed, registration.startScale + 1);
        registration.moving = create_multi_level_representation(...
            registration.moving, registration.startScale + 1);
        if strcmp(registration.method,'phase') && ...
                ~isempty(registration.fixedCertainty)
            registration.fixedCertainty = create_multi_level_representation(...
                registration.fixedCertainty, registration.startScale + 1);
        else
            registration.fixedCertainty = cell(1,registration.startScale + 1);
        end
        if strcmp(registration.method,'phase') && ...
                ~isempty(registration.movingCertainty)
            registration.movingCertainty = create_multi_level_representation(...
                registration.movingCertainty, registration.startScale + 1);
        else
            registration.movingCertainty = cell(1,registration.startScale + 1);
        end
    otherwise
        error('Unknow multi-level mode set');
end

%% Initialize variables
if registration.displayFinal
    fprintf('Initialize variables\n');
end
currentScale = registration.startScale;
iteration = 1;
endAlignmentLoop = false;
changeScale = false;

registration.displacementField = cell(registration.dims,1);
for k = 1 : registration.dims
    registration.displacementField{k} = zeros(size(registration.fixed{currentScale + 1}));
end
registration.certainty = zeros(size(registration.fixed{currentScale + 1}));
if isfield(registration,'initialTransformationMatrix')
    registration.transformationMatrix = registration.initialTransformationMatrix;
else
    registration.transformationMatrix = eye(registration.dims + 1);
end

if registration.logDomain && strcmp(registration.transformationModel,'non-rigid')
    registration.logDisplacementField = cell(registration.dims,1);
    for k = 1 : registration.dims
        registration.logDisplacementField{k} = zeros(size(registration.fixed{currentScale + 1}));
    end
end

%***********************************************************
% Evaluate registration error before starting
%***********************************************************

if registration.evaluate
    registration.deformed = registration.moving{currentScale + 1};
    
    % MSE
    registration.results.MSE(1) = sum_of_squared_differences(...
        registration.deformed,registration.fixed{currentScale + 1}) / ...
        numel(registration.deformed);
    
    % NCC
    registration.results.NCC(1) = normalized_cross_correlation(...
        registration.deformed,registration.fixed{currentScale + 1});
    
    if registration.dims == 3
        % Convert displacement cell to matrix
        displacementFieldMatrix = single(field_convert(registration.displacementField));
        displacementFieldMatrix = permute(displacementFieldMatrix, [2 3 4 1]);
        
        % harmonic
        harmonic = field_harmonic(displacementFieldMatrix,single([1 1 1]));
        registration.results.harmonicEnergy(1) = mean(harmonic(:));
        
        % jacobian
        jacobian = field_jacobian(displacementFieldMatrix,single([1 1 1]));
        registration.results.minJacobian(1) = min(jacobian(:));
        registration.results.maxJacobian(1) = max(jacobian(:));
        registration.results.numberOfVoxelsWithNegJacobian(1) = sum(jacobian(:) < 0);
    end
end

%% Run alignment loop
if registration.displayFinal
    fprintf('Starting registration using:\n');
    fprintf('Method:                    %s\n',registration.method);
    fprintf('TransformationModel:       %s\n',registration.transformationModel);
    fprintf('Scales:                    %i - %i\n',registration.stopScale,registration.startScale);
    fprintf('Iterations per scale:      %i\n',registration.iterationsPerScale);
    fprintf('Multi modal:               %i\n',registration.multiModal);
    fprintf('Symmetric:                 %i\n',registration.symmetric);
    if strcmp(registration.transformationModel,'non-rigid')
        fprintf('Fluid regularization:      %i\n',registration.fluidRegularization);
        fprintf('Elastic regularization:    %i\n',registration.elasticRegularization);
        fprintf('Log-domain:                %i\n',registration.logDomain);
        if registration.logDomain
            fprintf('BCH approximation:         %s\n',registration.BCHmode);
        else
            fprintf('Accumulation method:       %s\n',registration.accumulationMethod);
            fprintf('Apply certainty:           %i\n',registration.applyCertainty);
        end
    end
end
if registration.displayFinal
    fprintf('----------------------------------------\n');
end
while ~endAlignmentLoop
    
    if registration.displayFinal
        fprintf('Iteration %i of %i for scale %i\n',...
            iteration,registration.iterationsPerScale(currentScale + 1),currentScale);
    end
    
    %***********************************************************
    % Displacement estimation
    %***********************************************************
    
    if registration.displayFinal
        fprintf('Estimate update field\n');
    end
    if registration.symmetric
        [update] = estimate_symmetric_update_field(currentScale);
    else
        [update] = estimate_update_field(currentScale);
    end
    
    %***********************************************************
    % Non-rigid
    %***********************************************************
    
    if strcmp(registration.transformationModel,'non-rigid')
        
        %***********************************************************
        % Fluid regularization
        %***********************************************************
        
        if registration.fluidRegularization
            if registration.displayFinal
                fprintf('Regularization of the update field\n');
            end
            if registration.symmetric
                % Note when symmetric we are working in the log-domain
                % which means that the use of certainty is currently not
                % relevant
                update.displacementForward = field_regularization(...
                    update.displacementForward,...
                    [],...
                    false,...
                    registration.fluidRegularizationData(currentScale + 1));
                update.displacementBackward = field_regularization(...
                    update.displacementBackward,...
                    [],...
                    false,...
                    registration.fluidRegularizationData(currentScale + 1));
            else
                update.displacement = field_regularization(...
                    update.displacement,...
                    update.certainty,...
                    registration.applyCertainty,...
                    registration.fluidRegularizationData(currentScale + 1));
            end
        end
        
        %***********************************************************
        % Displacement field accumulation
        %***********************************************************
        
        if registration.displayFinal
            fprintf('Accumulate displacement field\n');
        end
        if registration.symmetric
            % Note when symmetric we are working in the log-domain which
            % means that the use of certainty is currently not relevant
            logDisplacementFieldForward = baker_campbell_hausdorff_formula(...
                registration.logDisplacementField,...
                update.displacementForward,...
                registration.BCHmode);
            logDisplacementFieldBackward = baker_campbell_hausdorff_formula(...
                backward(registration.logDisplacementField),...
                update.displacementBackward,...
                registration.BCHmode);
            for k = 1 : registration.dims
                registration.logDisplacementField{k} = 1/2 * ...
                    (logDisplacementFieldForward{k} - logDisplacementFieldBackward{k});
            end
        else
            if registration.logDomain
                % Note when working in the log-domain than the use of
                % certainty is currently not relevant
                registration.logDisplacementField = baker_campbell_hausdorff_formula(...
                    registration.logDisplacementField,update.displacement,...
                    registration.BCHmode);
            else
                if registration.applyCertainty
                    registration.displacementField = add_update_field_to_acc_field_with_cert(...
                        update.displacement, registration.displacementField,...
                        update.certainty, registration.certainty,...
                        registration.accumulationMethod, registration.interpolation);
                else
                    registration.displacementField = add_update_field_to_acc_field(...
                        update.displacement, registration.displacementField,...
                        registration.accumulationMethod, registration.interpolation);
                end
            end
        end
        
        %***********************************************************
        % Elastic regularization
        %***********************************************************
        
        if registration.elasticRegularization
            if registration.displayFinal
                fprintf('Regularization of the accumulated displacement field\n');
            end
            if registration.logDomain
                % Note when working in the log-domain than the use of
                % certainty is currently not relevant
                registration.logDisplacementField = field_regularization(...
                    registration.logDisplacementField,...
                    [],...
                    false,...
                    registration.elasticRegularizationData(currentScale + 1));
            else
                registration.displacementField = field_regularization(...
                    registration.displacementField,...
                    registration.certainty,...
                    registration.applyCertainty,...
                    registration.elasticRegularizationData(currentScale + 1));
            end
        end
        
        % Convert from velocity field to vector field
        if registration.logDomain
            registration.displacementField = field_exponentiation(...
                registration.logDisplacementField);
        end
        
        %***********************************************************
        % Translation or affine
        %***********************************************************
    else
        
        % If rigid transformation model remove scaling
        if strcmp(registration.transformationModel,'rigid')
            if registration.symmetric
                [U,S,V] = svd(update.transformationMatrixForward(1:end-1,1:end-1));
                update.transformationMatrixForward(1:end-1,1:end-1) = U*transpose(V);
                [U,S,V] = svd(update.transformationMatrixBackward(1:end-1,1:end-1));
                update.transformationMatrixBackward(1:end-1,1:end-1) = U*transpose(V);
            else
                [U,S,V] = svd(update.transformationMatrix(1:end-1,1:end-1));
                update.transformationMatrix(1:end-1,1:end-1) = U*transpose(V);
            end
        end
        %***********************************************************
        % Add update
        %***********************************************************
        if registration.symmetric
            registration.transformationMatrix = ...
                1/2*(update.transformationMatrixForward + inv(update.transformationMatrixBackward)) * ...
                registration.transformationMatrix;
        else
            registration.transformationMatrix = ...
                update.transformationMatrix*registration.transformationMatrix;
        end
        
        %***********************************************************
        % Compute current displacement
        %***********************************************************
        registration.displacementField = transformation2displacement(...
            registration.transformationMatrix, ...
            size(registration.moving{currentScale + 1}));
    end
    
    
    %***********************************************************
    % Display results of current iteration
    %***********************************************************
    
    if registration.display
        switch registration.dims
            case 2
                display_results_registration2d(...
                    registration.moving{1}, ...
                    registration.fixed{1}, ...
                    registration.displacementField, ...
                    registration.originalSize, ...
                    registration.gamma);
            case 3
                display_results_registration3d(...
                    registration.moving{1}, ...
                    registration.fixed{1}, ...
                    registration.displacementField, ...
                    registration.originalSize,...
                    registration.gamma);
            otherwise
        end
    end
    
    if registration.recordMovie
        switch registration.dims
            case 2
                moving = registration.moving{1};
                moving = uint8(255*moving/max(moving(:)));
                fixed = registration.fixed{1};
                fixed = uint8(255*fixed/max(fixed(:)));
                deformed = deformation(registration.moving{1}, ...
                    resampler(registration.displacementField, ...
                    registration.originalSize, 'relativeValues', true), ...
                    'linear');
                deformed = uint8(255*deformed/max(deformed(:)));
                comb = deformed;
                comb(:,:,2) = fixed;
                comb(:,:,3) = zeros(size(fixed));
                im = [repmat(moving, [1 1 3]) repmat(fixed, [1 1 3]); ...
                    repmat(deformed, [1 1 3]) comb];
                frame = im2frame(im);
                writeVideo(vidObj, frame);
            case 3
        end
        
    end
    
    %***********************************************************
    % Determine whether to change scale or not
    %***********************************************************
    
    if (registration.iterationsPerScale(currentScale + 1) <= iteration)
        if registration.displayFinal
            fprintf('----------------------------------------\n');
        end
        changeScale = true;
    end
    
    %***********************************************************
    % Evaluate registration error
    %***********************************************************
    
    if registration.evaluate
        registration.deformed = deformation(registration.moving{currentScale + 1},...
            registration.displacementField, registration.interpolation);
        
        % MSE
        registration.results.MSE(end + 1) = sum_of_squared_differences(...
            registration.deformed,registration.fixed{currentScale + 1}) / ...
            numel(registration.deformed);
        
        % NCC
        registration.results.NCC(end + 1) = normalized_cross_correlation(...
            registration.deformed,registration.fixed{currentScale + 1});
        
        if registration.dims == 3
            % Convert displacement cell to matrix
            displacementFieldMatrix = single(field_convert(registration.displacementField));
            displacementFieldMatrix = permute(displacementFieldMatrix, [2 3 4 1]);
            
            % Harmonic Energy
            harmonic = field_harmonic(displacementFieldMatrix,single([1 1 1]));
            registration.results.harmonicEnergy(end + 1) = mean(harmonic(:));
            
            % Jacobian
            jacobian = field_jacobian(displacementFieldMatrix,single([1 1 1]));
            registration.results.minJacobian(end + 1) = min(jacobian(:));
            registration.results.maxJacobian(end + 1) = max(jacobian(:));
            registration.results.numberOfVoxelsWithNegJacobian(end + 1) = sum(jacobian(:) < 0);
        end
    end
    
    %***********************************************************
    % Change scale if necessary update scale and iteration parameters and
    % resample displacementField, certainty and deformed model
    %***********************************************************
    
    if changeScale
        % Remember number of iterations
        registration.numberOfIterations(currentScale + 1) = iteration;
        
        currentScale = currentScale - 1;
        
        if currentScale < registration.stopScale
            endAlignmentLoop = true;
        else
            if (strcmp(registration.multiLevelMode,'resolution') && ...
                    strcmp(registration.transformationModel,'non-rigid'))
                if registration.logDomain
                    registration.logDisplacementField = resampler(...
                        registration.logDisplacementField,...
                        size(registration.moving{currentScale + 1}), ...
                        'interpolation', registration.interpolation, 'relativeValues', true);
                    registration.displacementField = field_exponentiation(...
                        registration.logDisplacementField);
                else
                    registration.displacementField = resampler(...
                        registration.displacementField,...
                        size(registration.moving{currentScale + 1}), ...
                        'interpolation', registration.interpolation, 'relativeValues', true);
                end
                
                registration.certainty = resampler(...
                    registration.certainty,...
                    size(registration.moving{currentScale + 1}),...
                    'interpolation', registration.interpolation);
            elseif (strcmp(registration.multiLevelMode,'resolution'))
                ratio = size(registration.moving{currentScale + 1}) ./ ...
                    size(registration.moving{currentScale + 2});
                registration.transformationMatrix(1:registration.dims,end) = ...
                    registration.transformationMatrix(1:registration.dims,end).*...
                    ratio(:);
                
                registration.displacementField = transformation2displacement(...
                    registration.transformationMatrix, ...
                    size(registration.moving{currentScale + 1}));
            end
            iteration = 1;
            changeScale = false;
        end
    else
        iteration = iteration + 1;
    end
end

%%
%***************************************************************
% Resample results to finest scale if needed and compute the
% deformation
%***************************************************************

if (registration.stopScale ~= 0)
    if registration.logDomain
        registration.logDisplacementField = resampler(...
            registration.logDisplacementField,...
            registration.originalSize, ...
            'interpolation', registration.interpolation, 'relativeValues', true);
        registration.displacementField = field_exponentiation(...
            registration.logDisplacementField);
    else
        registration.displacementField = resampler(registration.displacementField,...
            registration.originalSize,...
            'interpolation', registration.interpolation, 'relativeValues', true);
    end
    
    registration.certainty = resampler(registration.certainty, ...
        registration.originalSize,...
        'interpolation',registration.interpolation);
end

registration.deformed = deformation(registration.moving{1},...
    registration.displacementField,registration.interpolation);

if registration.symmetric && registration.logDomain
    registration.inverseDisplacementField = field_exponentiation_inverse(...
            registration.logDisplacementField);
end
%%
%***********************************************************
% Display end results
%***********************************************************

if registration.displayFinal
    switch registration.dims
        case 2
            display_results_registration2d(registration.moving{1}, ...
                registration.fixed{1}, ...
                registration.displacementField, ...
                registration.originalSize, ...
                registration.gamma);
        case 3
            display_results_registration3d(registration.moving{1}, ...
                registration.fixed{1}, ...
                registration.displacementField, ...
                registration.originalSize, ...
                registration.gamma);
        otherwise
    end
end

if registration.recordMovie
    close(vidObj);
end