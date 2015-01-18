function [G,h] = build_G_h_quadratic2d(A,delta_b,transformationModel)
% BUILD_G_H_QUADRATIC2D Builds the equation system G*p = h
% 
% [G,h] = build_G_h_quadratic2d(A,delta_b,transformationModel)
%
% INPUT ARGUMENTS
% A                     - A = [A22 A21]
%                             [A12 A11]
% delta_b               - b_fixed - b_moving
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

sz = size(delta_b);

A11 = vec(A(:,:,4));
A12 = 2*vec(A(:,:,3));
A21 = 2*vec(A(:,:,2));
A22 = vec(A(:,:,1));
B1 = vec(delta_b(:,:,2));
B2 = vec(delta_b(:,:,1));

switch transformationModel
    case 'translation'
        G = zeros(2,2);
        h = zeros(2,1);

        G(1,1) = sum(A11.^2 + A21.^2);
        G(1,2) = sum(A11.*A12 + A21.*A22);
        G(2,1) = G(1,2);
        G(2,2) = sum(A12.^2 + A22.^2);
        h(1) = sum(A11.*B1 + A21.*B2);
        h(2) = sum(A12.*B1 + A22.*B2);
    case {'rigid','affine'}
        [x,y] = meshgrid(1:sz(2),1:sz(1));
        x = x - sz(2)/2 - 0.5;
        y = y - sz(1)/2 - 0.5;
        x = x(:);
        y = y(:);
        
        G = [sum((A11.^2 + A21.^2).*x.^2),  sum((A11.^2 + A21.^2).*x.*y),   sum((A11.*A12 + A21.*A22).*x.^2),   sum((A11.*A12 + A21.*A22).*x.*y),   sum((A11.^2 + A21.^2).*x),      sum((A11.*A12 + A21.*A22).*x); ...
             0,                             sum((A11.^2 + A21.^2).*y.^2),   sum((A11.*A12 + A21.*A22).*x.*y),   sum((A11.*A12 + A21.*A22).*y.^2),   sum((A11.^2 + A21.^2).*y),      sum((A11.*A12 + A21.*A22).*y); ...
             0,                             0,                              sum((A12.^2 + A22.^2).*x.^2),       sum((A12.^2 + A22.^2).*x.*y),       sum((A11.*A12 + A21.*A22).*x),  sum((A12.^2 + A22.^2).*x); ...
             0,                             0,                              0,                                  sum((A12.^2 + A22.^2).*y.^2),       sum((A11.*A12 + A21.*A22).*y),  sum((A12.^2 + A22.^2).*y); ...
             0,                             0,                              0,                                  0,                                  sum(A11.^2 + A21.^2),           sum((A11.*A12 + A21.*A22)); ...
             0,                             0,                              0,                                  0,                                  0,                              sum(A12.^2 + A22.^2)];
        
        h = [sum((A11.*B1 + A21.*B2).*x); ...
             sum((A11.*B1 + A21.*B2).*y); ...
             sum((A12.*B1 + A22.*B2).*x); ...
             sum((A12.*B1 + A22.*B2).*y); ...
             sum((A11.*B1 + A21.*B2)); ...
             sum((A12.*B1 + A22.*B2))];
         
        G = G + transpose(G) - diag(diag(G));
end