function [a11, a12, a13, a22, a23, a33, b1, b2 b3] = ...
    estimate_displacement_based_upon_least_square_solution3d(solution, variates, ...
    sizeImage, ...
    dk, ck, ...
    filterDirection, ...
    numberOfFilterDirections, ...
    t11, t12, t13, t22, t23, t33)
% ESTIMATE_DISPLACEMENT_BASED_UPON_LEAST_SQUARE_SOLUTION3D Sets up the equation system to solve for estimating a displacement field
%
% [a11, a12, a13, a22, a23, a33, b1, b2 b3] = ...
%    estimate_displacement_based_upon_least_square_solution2d(solution, variates, ...
%                                                             sizeImage, ...
%                                                             dk, ck, ...
%                                                             filterDirection, ...
%                                                             numberOfFilterDirections, ...
%                                                             t11, t12, t13, t22, t23, t33)
%
% INPUT ARGUMENTS
% solution                      - Solution to use (only 3 is available in 3D)
% variates                      - Number of variates
% sizeImage                     - Image size
% dk                            - Phase-difference for k different filter
%                                 directions
% ck                            - Corresponding certainties for k different
%                                 filter directions
% filterDirection               - k filter directions
% numberOfFilterDirections      - Number of filter directions
% t11                           - Tensor element (1,1)
% t12                           - Tensor element (1,2)
% t13                           - Tensor element (1,3)
% t22                           - Tensor element (2,2)
% t23                           - Tensor element (2,3)
% t33                           - Tensor element (3,3)
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% a11 							- Components of local equation system A*d = b
% a12 							- Components of local equation system A*d = b
% a13 							- Components of local equation system A*d = b
% a22 							- Components of local equation system A*d = b
% a23 							- Components of local equation system A*d = b
% a33 							- Components of local equation system A*d = b
% b1 							- Components of local equation system A*d = b
% b2 							- Components of local equation system A*d = b
% b3 							- Components of local equation system A*d = b

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

a11 = zeros(sizeImage);
a12 = zeros(sizeImage);
a13 = zeros(sizeImage);
a22 = zeros(sizeImage);
a23 = zeros(sizeImage);
a33 = zeros(sizeImage);
b1  = zeros(sizeImage);
b2  = zeros(sizeImage);
b3  = zeros(sizeImage);

for k = 1 : variates
    switch solution
        case 0 %************************************************
            
        case 1 %************************************************
            
        case 2 %************************************************
            
        case 3 %************************************************
            tt11 = t11(:,:,:,k).^2 + t12(:,:,:,k).^2 + t13(:,:,:,k).^2;
            tt12 = t11(:,:,:,k).*t12(:,:,:,k) + t12(:,:,:,k).*t22(:,:,:,k) + t13(:,:,:,k).*t23(:,:,:,k);
            tt13 = t11(:,:,:,k).*t13(:,:,:,k) + t12(:,:,:,k).*t23(:,:,:,k) + t13(:,:,:,k).*t33(:,:,:,k);
            tt22 = t12(:,:,:,k).^2 +t22(:,:,:,k).^2 + t23(:,:,:,k).^2;
            tt23 = t12(:,:,:,k).*t13(:,:,:,k) + t22(:,:,:,k).*t23(:,:,:,k) + t23(:,:,:,k).*t33(:,:,:,k);
            tt33 = t13(:,:,:,k).^2 + t23(:,:,:,k).^2 + t33(:,:,:,k).^2;
            for l = 1 : numberOfFilterDirections
                a11 = a11 + ck(:,:,:,l,k).*tt11;
                a12 = a12 + ck(:,:,:,l,k).*tt12;
                a13 = a13 + ck(:,:,:,l,k).*tt13;
                a22 = a22 + ck(:,:,:,l,k).*tt22;
                a23 = a23 + ck(:,:,:,l,k).*tt23;
                a33 = a33 + ck(:,:,:,l,k).*tt33;
                
                b1  =  b1 + ck(:,:,:,l,k).*dk(:,:,:,l,k).*(filterDirection{l}(1)*tt11 + (filterDirection{l}(2)*tt12) + (filterDirection{l}(3)*tt13));
                b2  =  b2 + ck(:,:,:,l,k).*dk(:,:,:,l,k).*(filterDirection{l}(1)*tt12 + (filterDirection{l}(2)*tt22) + (filterDirection{l}(3)*tt23));
                b3  =  b3 + ck(:,:,:,l,k).*dk(:,:,:,l,k).*(filterDirection{l}(1)*tt13 + (filterDirection{l}(2)*tt23) + (filterDirection{l}(3)*tt33));
            end
        case 4 %************************************************
            
    end
end