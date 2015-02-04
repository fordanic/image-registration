function [updateField, updateCertainty] = ...
    morphon3d(moving, movingCertainty, fixed, fixedCertainty, ...
    scale, varargin)
% MORPHON3D Estimates a displacement field between moving and fixed using the Morphon
%
% [updateField, updateCertainty] = morphon3d(...
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
% 'tensorSigma'         - Sigma of Gaussian to smooth tensor (2.5)
% 'elementsSigma'       - Sigma of Gaussian to smooth equation system (2.5)
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
tensorSigma = 2.5;
elementsSigma = 2.5;

for k=1:2:length(varargin),         % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

%***************************************************************
% Load quadrature filters
%***************************************************************

load quadratureFiltersForMorphonRegistration3D

numberOfFilterDirections = 6;

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
t13 = zeros(sizeImage);
t22 = zeros(sizeImage);
t23 = zeros(sizeImage);
t33 = zeros(sizeImage);
updateField = cell(3,1);

%***************************************************************
% Compute ck, dk and T
%***************************************************************

for k = 1 : numberOfFilterDirections
    % Filter moving and fixed
    qa = imfilter(moving, f{k}, 'same', 'conv');
    qb = imfilter(fixed, f{k}, 'same', 'conv');
    
    % Estimate dk and ck
    qq = qa.*conj(qb);
    Aqq = abs(qq);
    dk(:,:,:,k) = atan2(imag(qq),real(qq));
    ck(:,:,:,k) = Aqq;
    ck(:,:,:,k) = ck(:,:,:,k).^0.5.*(cos(dk(:,:,:,k)/2).^2);
    
    % Estimate T (Ta and Tb)
    Aqb = abs(qb);
    t11 = t11 + Aqb*m11{k};
    t12 = t12 + Aqb*m12{k};
    t13 = t13 + Aqb*m13{k};
    t22 = t22 + Aqb*m22{k};
    t23 = t23 + Aqb*m23{k};
    t33 = t33 + Aqb*m33{k};
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
    Tcert = sqrt(t11.^2 + 2*t12.^2 + 2*t13.^2 + t22.^2 + 2*t23.^2 + t33.^2);
    t11 = normgauss_smoothing(t11, Tcert, tensorSigma);
    t12 = normgauss_smoothing(t12, Tcert, tensorSigma);
    t13 = normgauss_smoothing(t13, Tcert, tensorSigma);
    t22 = normgauss_smoothing(t22, Tcert, tensorSigma);
    t23 = normgauss_smoothing(t23, Tcert, tensorSigma);
    t33 = normgauss_smoothing(t33, Tcert, tensorSigma);
else
    t11 = gauss_smoothing(t11, tensorSigma);
    t12 = gauss_smoothing(t12, tensorSigma);
    t13 = gauss_smoothing(t13, tensorSigma);
    t22 = gauss_smoothing(t22, tensorSigma);
    t23 = gauss_smoothing(t23, tensorSigma);
    t33 = gauss_smoothing(t33, tensorSigma);
end

tNorm = (sqrt(t11.^2 + 2*t12.^2 + 2*t13.^2 + t22.^2 + 2*t23.^2 + t33.^2) + eps);
t11 = t11./max(tNorm(:));
t12 = t12./max(tNorm(:));
t13 = t13./max(tNorm(:));
t22 = t22./max(tNorm(:));
t23 = t23./max(tNorm(:));
t33 = t33./max(tNorm(:));

%***************************************************************
% Solutions based upon minimizing a least square problem, 
% i.e. runParameters.solution = 0-5
%***************************************************************

if (solution >= 0 && solution <= 5)
    
    [a11, a12, a13, a22, a23, a33, b1, b2, b3] = ...
        estimate_displacement_based_upon_least_square_solution3d(...
            solution, 1, ...
            sizeImage, ...
            dk, ck, ...
            filterDirection, ...
            numberOfFilterDirections, ...
            t11, t12, t13, t22, t23, t33);

    %***********************************************************
    % Reqularize axx and bx
    %***********************************************************
    
    a11 = gauss_smoothing(a11, elementsSigma);
    a12 = gauss_smoothing(a12, elementsSigma);
    a13 = gauss_smoothing(a13, elementsSigma);
    a22 = gauss_smoothing(a22, elementsSigma);
    a23 = gauss_smoothing(a23, elementsSigma);
    a33 = gauss_smoothing(a33, elementsSigma);
    
    b1 = gauss_smoothing(b1, elementsSigma);
    b2 = gauss_smoothing(b2, elementsSigma);
    b3 = gauss_smoothing(b3, elementsSigma);
    
    %***********************************************************
    % Solve equation system and compute incremental displacement and 
    % incremental certainty
    %***********************************************************
    
    updateCertainty = a11 + a22 + a33; 

    norm = 1./(a11.*a22.*a33 - a11.*a23.*a23 - a12.*a12.*a33 + a12.*a23.*a13 + a13.*a12.*a23 - a13.*a22.*a13 + eps);
    
    updateField{2} = norm.*((b3.*(a12.*a23 - a13.*a22)) - (b2.*(a12.*a33 - a13.*a23)) + (b1.*(a22.*a33 - a23.*a23)));
    updateField{1} = norm.*((b2.*(a11.*a33 - a13.*a13)) - (b3.*(a11.*a23 - a13.*a12)) - (b1.*(a12.*a33 - a23.*a13)));
    updateField{3} = norm.*((b3.*(a11.*a22 - a12.*a12)) - (b2.*(a11.*a23 - a12.*a13)) + (b1.*(a12.*a23 - a22.*a13)));
    
end

%***************************************************************
% Step size compensation
%***************************************************************

% Scale independent ampfactor
displacementAmplification = 2;

for k = 1 : 3
    updateField{k} = ...
        -displacementAmplification*updateField{k};
end

% Fix to gain more power on fine detail scales
updateCertainty = updateCertainty.*2^-scale;  
