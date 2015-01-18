function [G,h] = build_G_h_quadratic3d(A,delta_b,transformationModel)
% BUILD_G_H_QUADRATIC3D Builds the equation system G*p = h
% 
%
% INPUT ARGUMENTS
% b1                    - b = [b1 b2]'
% b2                    - b = [b1 b2]'
% deltaC                - c_fixed - c_moving
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

sz = size(delta_b(:,:,:,1));

A11 = vec(A(:,:,:,5));
A12 = 2*vec(A(:,:,:,4));
A13 = 2*vec(A(:,:,:,8));
A21 = 2*vec(A(:,:,:,2));
A22 = vec(A(:,:,:,1));
A23 = 2*vec(A(:,:,:,7));
A31 = 2*vec(A(:,:,:,6));
A32 = 2*vec(A(:,:,:,3));
A33 = vec(A(:,:,:,9));
B1 = vec(delta_b(:,:,:,2));
B2 = vec(delta_b(:,:,:,1));
B3 = vec(delta_b(:,:,:,3));

switch transformationModel
    case 'translation'
        G = zeros(3,3);
        h = zeros(3,1);

        G(1,1) = sum(A11.^2 + A21.^2 + A31.^2);
        G(1,2) = sum(A11.*A12 + A21.*A22 + A31.*A32);
        G(1,3) = sum(A11.*A13 + A21.*A23 + A31.*A33);
        G(2,1) = G(1,2);
        G(2,2) = sum(A12.^2 + A22.^2 + A32.^2);
        G(2,3) = sum(A12.*A13 + A22.*A23 + A32.*A33);
        G(3,1) = G(1,3);
        G(3,2) = G(2,3);
        G(3,3) = sum(A13.^2 + A23.^2 + A33.^2);
        h(1) = sum(A11.*B1 + A21.*B2 + A31.*B3);
        h(2) = sum(A12.*B1 + A22.*B2 + A32.*B3);
        h(3) = sum(A13.*B1 + A23.*B2 + A33.*B3);
    case {'rigid','affine'}
        [x,y,z] = meshgrid(1:sz(2),1:sz(1),1:sz(3));
        x = x - sz(2)/2 - 0.5;
        y = y - sz(1)/2 - 0.5;
        z = z - sz(3)/2 - 0.5;
        x = x(:);
        y = y(:);
        z = z(:);
        G = [   sum((A11.^2 + A21.^2 + A31.^2).*x.^2),  sum((A11.^2 + A21.^2 + A31.^2).*x.*y),  sum((A11.^2 + A21.^2 + A31.^2).*x.*z),  sum((A11.*A12 + A21.*A22 + A31.*A32).*x.^2),    sum((A11.*A12 + A21.*A22 + A31.*A32).*x.*y),    sum((A11.*A12 + A21.*A22 + A31.*A32).*x.*z),    sum((A11.*A13 + A21.*A23 + A31.*A33).*x.^2),    sum((A11.*A13 + A21.*A23 + A31.*A33).*x.*y),    sum((A11.*A13 + A21.*A23 + A31.*A33).*x.*z),    sum((A11.^2 + A21.^2 + A31.^2).*x),         sum((A11.*A12 + A21.*A22 + A31.*A32).*x),   sum((A11.*A13 + A21.*A23 + A31.*A33).*x); ...
                0,                                      sum((A11.^2 + A21.^2 + A31.^2).*y.^2),  sum((A11.^2 + A21.^2 + A31.^2).*y.*z),  sum((A11.*A12 + A21.*A22 + A31.*A32).*x.*y),    sum((A11.*A12 + A21.*A22 + A31.*A32).*y.^2),    sum((A11.*A12 + A21.*A22 + A31.*A32).*y.*z),    sum((A11.*A13 + A21.*A23 + A31.*A33).*x.*y),    sum((A11.*A13 + A21.*A23 + A31.*A33).*y.^2),    sum((A11.*A13 + A21.*A23 + A31.*A33).*y.*z),    sum((A11.^2 + A21.^2 + A31.^2).*y),         sum((A11.*A12 + A21.*A22 + A31.*A32).*y),   sum((A11.*A13 + A21.*A23 + A31.*A33).*y); ...
                0,                                      0,                                      sum((A11.^2 + A21.^2 + A31.^2).*z.^2),  sum((A11.*A12 + A21.*A22 + A31.*A32).*x.*z),    sum((A11.*A12 + A21.*A22 + A31.*A32).*y.*z),    sum((A11.*A12 + A21.*A22 + A31.*A32).*z.^2),    sum((A11.*A13 + A21.*A23 + A31.*A33).*x.*z),    sum((A11.*A13 + A21.*A23 + A31.*A33).*y.*z),    sum((A11.*A13 + A21.*A23 + A31.*A33).*z.^2),    sum((A11.^2 + A21.^2 + A31.^2).*z),         sum((A11.*A12 + A21.*A22 + A31.*A32).*z),   sum((A11.*A13 + A21.*A23 + A31.*A33).*z); ...
                0,                                      0,                                      0,                                      sum((A12.^2 + A22.^2 + A32.^2).*x.^2),          sum((A12.^2 + A22.^2 + A32.^2).*x.*y),          sum((A12.^2 + A22.^2 + A32.^2).*x.*z),          sum((A12.*A13 + A22.*A23 + A32.*A33).*x.^2),    sum((A12.*A13 + A22.*A23 + A32.*A33).*x.*y),    sum((A12.*A13 + A22.*A23 + A32.*A33).*x.*z),    sum((A11.*A12 + A21.*A22 + A31.*A32).*x),   sum((A12.^2 + A22.^2 + A32.^2).*x),         sum((A12.*A13 + A22.*A23 + A32.*A33).*x); ...
                0,                                      0,                                      0,                                      0,                                              sum((A12.^2 + A22.^2 + A32.^2).*y.^2),          sum((A12.^2 + A22.^2 + A32.^2).*y.*z),          sum((A12.*A13 + A22.*A23 + A32.*A33).*x.*y),    sum((A12.*A13 + A22.*A23 + A32.*A33).*y.^2),    sum((A12.*A13 + A22.*A23 + A32.*A33).*y.*z),    sum((A11.*A12 + A21.*A22 + A31.*A32).*y),   sum((A12.^2 + A22.^2 + A32.^2).*y),         sum((A12.*A13 + A22.*A23 + A32.*A33).*y); ...
                0,                                      0,                                      0,                                      0,                                              0,                                              sum((A12.^2 + A22.^2 + A32.^2).*z.^2),          sum((A12.*A13 + A22.*A23 + A32.*A33).*x.*z),    sum((A12.*A13 + A22.*A23 + A32.*A33).*y.*z),    sum((A12.*A13 + A22.*A23 + A32.*A33).*z.^2),    sum((A11.*A12 + A21.*A22 + A31.*A32).*z),   sum((A12.^2 + A22.^2 + A32.^2).*z),         sum((A12.*A13 + A22.*A23 + A32.*A33).*z); ...
                0,                                      0,                                      0,                                      0,                                              0,                                              0,                                              sum((A13.^2 + A23.^2 + A33.^2).*x.^2),          sum((A13.^2 + A23.^2 + A33.^2).*x.*y),          sum((A13.^2 + A23.^2 + A33.^2).*x.*z),          sum((A11.*A13 + A21.*A23 + A31.*A33).*x),   sum((A12.*A13 + A22.*A23 + A32.*A33).*x),   sum((A13.^2 + A23.^2 + A33.^2).*x); ...
                0,                                      0,                                      0,                                      0,                                              0,                                              0,                                              0,                                              sum((A13.^2 + A23.^2 + A33.^2).*y.^2),          sum((A13.^2 + A23.^2 + A33.^2).*y.*z),          sum((A11.*A13 + A21.*A23 + A31.*A33).*y),   sum((A12.*A13 + A22.*A23 + A32.*A33).*y),   sum((A13.^2 + A23.^2 + A33.^2).*y); ...
                0,                                      0,                                      0,                                      0,                                              0,                                              0,                                              0,                                              0,                                              sum((A13.^2 + A23.^2 + A33.^2).*z.^2),          sum((A11.*A13 + A21.*A23 + A31.*A33).*z),   sum((A12.*A13 + A22.*A23 + A32.*A33).*z),   sum((A13.^2 + A23.^2 + A33.^2).*z); ...
                0,                                      0,                                      0,                                      0,                                              0,                                              0,                                              0,                                              0,                                              0,                                              sum((A11.^2 + A21.^2 + A31.^2)),            sum((A11.*A12 + A21.*A22 + A31.*A32)),      sum((A11.*A13 + A21.*A23 + A31.*A33)); ...
                0,                                      0,                                      0,                                      0,                                              0,                                              0,                                              0,                                              0,                                              0,                                              0,                                          sum((A12.^2 + A22.^2 + A32.^2)),            sum((A12.*A13 + A22.*A23 + A32.*A33)); ...
                0,                                      0,                                      0,                                      0,                                              0,                                              0,                                              0,                                              0,                                              0,                                              0,                                          0,                                          sum((A13.^2 + A23.^2 + A33.^2))];
 
        h = [sum((A11.*B1 + A21.*B2 + A31.*B3).*x); ...
             sum((A11.*B1 + A21.*B2 + A31.*B3).*y); ...
             sum((A11.*B1 + A21.*B2 + A31.*B3).*z); ...
             sum((A12.*B1 + A22.*B2 + A32.*B3).*x); ...
             sum((A12.*B1 + A22.*B2 + A32.*B3).*y); ...
             sum((A12.*B1 + A22.*B2 + A32.*B3).*z); ...
             sum((A13.*B1 + A23.*B2 + A33.*B3).*x); ...
             sum((A13.*B1 + A23.*B2 + A33.*B3).*y); ...
             sum((A13.*B1 + A23.*B2 + A33.*B3).*z); ...
             sum((A11.*B1 + A21.*B2 + A31.*B3)); ...
             sum((A12.*B1 + A22.*B2 + A32.*B3)); ...
             sum((A13.*B1 + A23.*B2 + A33.*B3))];
         
        G = G + transpose(G) - diag(diag(G));
end
