function [A,h] = build_A_h_2d(dphi,cert,phiGrad,transformatioinModel)
% BUILD_A_H_LINEAR2D Builds the equation system A*p = h
%
% INPUT ARGUMENTS
% dphi                  - Phase-difference
% cert                  - Certainty
% phiGrad               - Phase gradient
% transformationModel   - Transformation model (translation or affine)
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% A 					- A matrix
% h 					- h vector
%
% See "Phase-Based Multidimensional Volume Registration" by Hemmendorf et al
% or "PHASE BASED VOLUME REGISTRATION USING CUDA" by Eklund et al for a
% detailed description of how the equation system is set up. Note that here
% we use a differenct S(x)p then the one used in the papers. This is done
% in order to be consistent with functions build_G_h_linear2d and
% build_G_h_linear3d.

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

sz = size(dphi(:,:,1));

dphiX = vec(dphi(:,:,1));
dphiY = vec(dphi(:,:,2));
certX = vec(cert(:,:,1));
certY = vec(cert(:,:,2));
phiGradX = vec(phiGrad(:,:,1));
phiGradY = vec(phiGrad(:,:,2));

switch transformatioinModel
    case 'translation'
        A = [   sum(certX.*phiGradX.^2),    0; ...
                0,                          sum(certY.*phiGradY.^2)];
        h = [   sum(certX.*dphiX.*phiGradX); ...
                sum(certY.*dphiY.*phiGradY)];
    case {'rigid','affine'}
        [x,y] = meshgrid(1:sz(2),1:sz(1));
        x = x - sz(2)/2 - 0.5;
        y = y - sz(1)/2 - 0.5;
        x = x(:);
        y = y(:);
        
        A = [   sum(certX.*phiGradX.^2.*x.^2),  sum(certX.*phiGradX.^2.*x.*y),  0,                                  0,                              sum(certX.*phiGradX.^2.*x), 0; ...
                0,                              sum(certX.*phiGradX.^2.*y.^2),  0,                                  0,                              sum(certX.*phiGradX.^2.*y), 0; ...
                0,                              0,                              sum(certY.*phiGradY.^2.*x.^2),      sum(certY.*phiGradY.^2.*x.*y),  0,                          sum(certY.*phiGradY.^2.*x); ...
                0,                              0,                              0,                                  sum(certY.*phiGradY.^2.*y.^2),  0,                          sum(certY.*phiGradY.^2.*y); ...
                0,                              0,                              0,                                  0,                              sum(certX.*phiGradX.^2),    0; ...
                0,                              0,                              0,                                  0,                              0,                          sum(certY.*phiGradY.^2)];
        
        h = [   sum(certX.*dphiX.*phiGradX.*x); ...
                sum(certX.*dphiX.*phiGradX.*y); ...
                sum(certY.*dphiY.*phiGradY.*x); ...
                sum(certY.*dphiY.*phiGradY.*y); ...
                sum(certX.*dphiX.*phiGradX); ...
                sum(certY.*dphiY.*phiGradY)];
        
        A = A + transpose(A) - diag(diag(A));
end