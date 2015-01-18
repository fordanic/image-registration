function [G,h] = build_G_h_linear2d(b,delta_c,transformationModel)
% BUILD_G_H_LINEAR2D Builds the equation system G*p = h
% 
% INPUT ARGUMENTS
% b                     - b = [b2 b1]'
% delta_c               - c_fixed - c_moving
% transformationModel   - Transformation model (translation or affine)
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% G 					- G matrix
% h 					- h vector
%
% See "Affine and Deformable Registration Based on Polynomial Expansion"
% by Gunnar Farnebäck and Carl-Fredrik Westin for detailed exaplanation on
% how the equation system is set up

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

sz = size(b);

b1 = vec(b(:,:,2));
b2 = vec(b(:,:,1));
delta_c = vec(delta_c);

switch transformationModel
    case 'translation'
        G = zeros(2,2);
        h = zeros(2,1);

        G(1,1) = sum(b1.*b1);
        G(1,2) = sum(b1.*b2);
        G(2,1) = G(1,2);
        G(2,2) = sum(b2.*b2);
        h(1) = sum(b1.*delta_c);
        h(2) = sum(b2.*delta_c);
    case {'rigid','affine'}
        [x,y] = meshgrid(1:sz(2),1:sz(1));
        x = x - sz(2)/2 - 0.5;
        y = y - sz(1)/2 - 0.5;
        x = x(:);
        y = y(:);
        
        G = [sum(b1.^2.*x.^2),  sum(b1.^2.*x.*y),   sum(b1.*b2.*x.^2),  sum(b1.*b2.*x.*y),  sum(b1.^2.*x),  sum(b1.*b2.*x); ...
             0,                 sum(b1.^2.*y.^2),   sum(b1.*b2.*x.*y),  sum(b1.*b2.*y.^2),  sum(b1.^2.*y),  sum(b1.*b2.*y); ...
             0,                 0,                  sum(b2.^2.*x.^2),   sum(b2.^2.*x.*y),   sum(b1.*b2.*x), sum(b2.^2.*x); ...
             0,                 0,                  0,                  sum(b2.^2.*y.^2),   sum(b1.*b2.*y), sum(b2.^2.*y); ...
             0,                 0,                  0,                  0,                  sum(b1.^2),     sum(b1.*b2); ...
             0,                 0,                  0,                  0,                  0,              sum(b2.^2)];
        
        h = [sum(b1.*delta_c.*x); ...
             sum(b1.*delta_c.*y); ...
             sum(b2.*delta_c.*x); ...
             sum(b2.*delta_c.*y); ...
             sum(b1.*delta_c); ...
             sum(b2.*delta_c)];
         
        G = G + transpose(G) - diag(diag(G));
end