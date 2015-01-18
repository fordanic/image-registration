function check_registration()
% CHECK_REGISTRATION Checks that the data struct registration is ok
%
% check_registration
%
% INPUT ARGUMENTS
% N/A
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% N/A

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

global registration;

% Check data size
dims = ndims(registration.fixed);

sizeFixed = size(registration.fixed);
sizeMoving = size(registration.moving);

if sum(sizeFixed == sizeMoving) ~= dims
    error('Provided image data is not of the same size')
end

% Set to deafult false
registration.useTensorMagnitudeRepresentation = false;

if ~isfield(registration,'multiModal')
    registration.multiModal = false;
end
if ~isfield(registration,'symmetric')
    registration.symmetric = false;
end
if ~isfield(registration,'logDomain')
    registration.logDomain = false;
end

% Check multi-level mode
if ~isfield(registration,'multiLevelMode')
    registration.multiLevelMode = 'resolution';
end
if ~strcmp(registration.multiLevelMode,'scale') && ~strcmp(registration.multiLevelMode,'resolution')
    error('Unknown multi level mode');
end

% Check interpolation
if ~isfield(registration,'interpolation')
    registration.interpolation = 'linear';
end

% Check scale parameters
if ~isfield(registration,'startScale')
    error('A start scale must be provided')
end
if ~isfield(registration,'stopScale')
    error('A stop scale must be provided')
end
if registration.startScale < registration.stopScale
    error('Stop scale must be smaller than start scale')
end
if ~isfield(registration,'iterationsPerScale')
    error('Number of iterations per scale must be defined')
end
if (length(registration.iterationsPerScale) == 1)
    registration.iterationsPerScale = ...
        registration.iterationsPerScale * ones(1,registration.startScale + 1);
elseif (length(registration.iterationsPerScale) ~= registration.startScale + 1)
    error('Mismatch for data provided for min iteration per scale');
end

%% Check method specific parameters
switch registration.method
    case 'phase'
        if registration.multiModal
            registration.useTensorMagnitudeRepresentation = true;
        end
        if isfield(registration,'fixedCertainty') && ...
                ~isempty(registration.fixedCertainty)
            sizeFixedCertainty = size(registration.fixedCertainty);
            if sum(sizeFixed == sizeFixedCertainty) ~= dims
                error('Provided certainty mask of fixed image is not of same size as fixed image')
            end
        else
            registration.fixedCertainty = [];
        end
        if isfield(registration,'movingCertainty') && ...
                ~isempty(registration.movingCertainty)
            sizeMovingCertainty = size(registration.movingCertainty);
            if sum(sizeMoving == sizeMovingCertainty) ~= dims
                error('Provided certainty mask of moving image is not of same size as moving image')
            end
        else
            registration.movingCertainty = [];
        end
    case 'optical-flow'
        if registration.multiModal
            if (length(registration.numberOfChannels) == 1)
                registration.numberOfChannels = ...
                    registration.numberOfChannels * ones(1,registration.startScale + 1);
            elseif (length(registration.numberOfChannels) ~= registration.startScale + 1)
                error('Mismatch for data provided for min iteration per scale');
            end
        end
    case 'polynomial-expansion'
        if ~isfield(registration,'signalModel')
            warning('Signal model is not set, setting to linear');
            registration.signalModel = 'linear';
        end
        if registration.multiModal
            if (length(registration.numberOfChannels) == 1)
                registration.numberOfChannels = ...
                    registration.numberOfChannels * ones(1,registration.startScale + 1);
            elseif (length(registration.numberOfChannels) ~= registration.startScale + 1)
                error('Mismatch for data provided for min iteration per scale');
            end
            
            if strcmp(registration.signalModel,'quadratic')
                warning('Multi-modal is not supported for a quadratic signal model in polynomial expansion, setting to false')
                registration.multiModal = false;
            end
        end
    otherwise
        error('Unknown registration method')
end

%% Check transformation model parameters
if strcmp(registration.transformationModel,'non-rigid')
    % Check regularization parameters
    
    % Fluid regularization
    if registration.fluidRegularization
        if (length(registration.fluidRegularizationData) == 1)
            registration.fluidRegularizationData = ...
                registration.fluidRegularizationData * ones(1,registration.startScale + 1);
        elseif (length(registration.fluidRegularizationData) ~= registration.startScale + 1)
            error('Mismatch for data provided for fluid regularization');
        end
    end
    
    % Elastic regularization
    if registration.elasticRegularization
        if (length(registration.elasticRegularizationData) == 1)
            registration.elasticRegularizationData = ...
                registration.elasticRegularizationData * ones(1,registration.startScale + 1);
        elseif (length(registration.elasticRegularizationData) ~= registration.startScale + 1)
            error('Mismatch for data provided for elastic regularization');
        end
    end
    
    if registration.symmetric
        registration.logDomain = true;
    end
    
    if ~isfield(registration,'applyCertainty')
        registration.applyCertainty = false;
    end
    
    if registration.logDomain
        % Check BCHmode
        if ~isfield(registration,'BCHmode')
            warning('BCH approximation mode is note set, setting to A');
            registration.BCHmode = 'A';
        end
        % Check apply certainty
        if registration.applyCertainty
            warning('The use of certainty is not relevant when working in the log-domain, setting to false')
            registration.applyCertainty = false;
        end
    end
end

if ~isfield(registration,'evaluate')
    registration.evaluate = false;
end

if ~isfield(registration,'display')
    registration.display = false;
end

if ~isfield(registration,'displayFinal')
    registration.displayFinal = false;
end

if ~isfield(registration,'gamma')
    registration.gamma = 0.5;
end

if ~isfield(registration,'recordMovie')
    registration.recordMovie = false;
end

if dims == 3 && registration.recordMovie == true
    warning('Movie recording is only supported for 2D data, setting to false')
    registration.recordMovie = false;
end