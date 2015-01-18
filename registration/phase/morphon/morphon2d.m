function [updateField, updateCertainty] = ...
    morphon2d(moving, movingCertainty, fixed, fixedCertainty, ...
    scale, varargin)
% MORPHON2D Estimates a displacement field between moving and fixed using the Morphon
%
% [updateField, updateCertainty] = morphon2d(...
%     moving, movingCertainty, fixed, fixedCertainty, scale)
%
% INPUT ARGUMENTS
% moving                - Moving image
% movingCertainty       - Certainty mask of moving image
% fixed                 - Fixed image
% fixedCertainty        - Certainty mask of fixed image
% scale                 - Current scale
%
% OPTIONAL INPUT ARGUMENTS
% 'solution'            - Solution to apply when estimating the
%                         displacement field (3)
% 'applyCertainty'      - Apply certainty when smoothing tensor elements
%                         and equation system elements (true)
% 'tensorSigma'         - Sigma of Gaussian to smooth tensor (1.5)
% 'elementsSigma'       - Sigma of Gaussian to smooth equation system (1.5)
%
% OUTPUT ARGUMENTS
% updateField 			- Estimated displacement
% updateCertainty 		- Certainty of displacement

% Copyright (c) 2011 Daniel Forsberg
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

% Setup default parameters
solution = 3;
applyCertainty = true;
tensorSigma = 1.5;
elementsSigma = 1.5;

% Overwrites default parameter
for k=1:2:length(varargin)
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

%***************************************************************
% Load quadrature filters
%***************************************************************

load quadratureFiltersForMorphonRegistration2D

numberOfFilterDirections = 4;

%***************************************************************
% Determine image dimension
%***************************************************************

sizeImage = size(fixed);

%***************************************************************
% Initialize some variables
%***************************************************************

dk = zeros([sizeImage numberOfFilterDirections]);
ck = zeros([sizeImage numberOfFilterDirections]);
t11 = zeros(sizeImage);
t12 = zeros(sizeImage);
t22 = zeros(sizeImage);
updateField = cell(2,1);

%***************************************************************
% Compute ck, dk and T
%***************************************************************

for k = 1:numberOfFilterDirections
    % Filter moving and fixed
    qa = imfilter(moving, f{k}, 'same', 'conv');
    qb = imfilter(fixed, f{k}, 'same', 'conv');
    
    % Estimate dk and ck
    qq = qa.*conj(qb);
    Aqq = abs(qq);
    dk(:,:,k) = atan2(imag(qq),real(qq));
    ck(:,:,k) = Aqq;
    ck(:,:,k) = ck(:,:,k).^0.5.*(cos(dk(:,:,k)/2).^2);
    
    % Estimate T (Ta and Tb)
    Aqb = abs(qb);
    t11 = t11 + Aqb*m11{k};
    t12 = t12 + Aqb*m12{k};
    t22 = t22 + Aqb*m22{k};
end

% Add certainty masks
if ~isempty(fixedCertainty)
    ck = bsxfun(@times,ck,fixedCertainty);
end
if ~isempty(movingCertainty)
    ck = bsxfun(@times,ck,movingCertainty);
end

%***************************************************************
% Average and normalize Ta and Tb
%***************************************************************

if applyCertainty
    Tcert = sqrt(t11.^2 +t22.^2 + 2*t12.^2);
    t11 = normgauss_smoothing(t11, Tcert, tensorSigma);
    t12 = normgauss_smoothing(t12, Tcert, tensorSigma);
    t22 = normgauss_smoothing(t22, Tcert, tensorSigma);
else
    t11 = gauss_smoothing(t11, tensorSigma);
    t12 = gauss_smoothing(t12, tensorSigma);
    t22 = gauss_smoothing(t22, tensorSigma);
end

tNorm = (sqrt(t11.^2 +t22.^2 + 2*t12.^2) + eps);
t11 = t11./max(tNorm(:));
t12 = t12./max(tNorm(:));
t22 = t22./max(tNorm(:));

%***************************************************************
% Solutions based upon minimizing a least square problem, 
% i.e. runParameters.solution = 0-5
%***************************************************************

if (solution >= 0 && solution <= 5)
    
    [a11, a12, a21, a22, b1, b2] = ...
        estimate_displacement_based_upon_least_square_solution2d(...
            solution, ...
            sizeImage, ...
            dk, ck, ...
            filterDirection, ...
            numberOfFilterDirections, ...
            t11, t12, t22);

    %***********************************************************
    % Reqularize axx and bx
    %***********************************************************
    
    a11 = gauss_smoothing(a11, elementsSigma);
    a12 = gauss_smoothing(a12, elementsSigma);
    a21 = gauss_smoothing(a21, elementsSigma);
    a22 = gauss_smoothing(a22, elementsSigma);
    
    b1 = gauss_smoothing(b1, elementsSigma);
    b2 = gauss_smoothing(b2, elementsSigma);
    
    %***********************************************************
    % Solve equation system and compute incremental displacement and 
    % incremental certainty
    %***********************************************************
    
    % Determine eigenvalues of matrix A
    root = sqrt((a11-a22).^2+4.*a12.^2);
    if ~isreal(root)
        error('The eigenvalues are complex!')
    end
    
    %lambda1 = (a11+a22+root)./2;
    lambda2 = abs((a11+a22-root)./2);
    
    updateCertainty = lambda2;
    
    %updateCertainty = a11 + a22; % use trace as certainty

    norm = 1./(a11.*a22 - a12.*a21 + eps);
    
    %updateCertainty = 1./norm;
    
    updateField{2} = norm.*( a22.*b1 - a12.*b2);
    updateField{1} = norm.*(-a21.*b1 + a11.*b2);
    
end

%***************************************************************
% Step size compensation
%***************************************************************

% Scale independent ampfactor
displacementAmplification = 2;

for k = 1 : 2
    updateField{k} = ...
        -displacementAmplification*updateField{k};
end

% Fix to gain more power on fine detail scales
updateCertainty = updateCertainty.*2^-scale;
