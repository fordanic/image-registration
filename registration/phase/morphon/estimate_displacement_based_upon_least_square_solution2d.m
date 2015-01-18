function [a11, a12, a21, a22, b1, b2] = ...
    estimate_displacement_based_upon_least_square_solution2d(solution, ...
                                                             sizeImage, ...
                                                             dk, ck, ...
                                                             filterDirection, ...
                                                             numberOfFilterDirections, ...
                                                             t11, t12, t22)
% ESTIMATE_DISPLACEMENT_BASED_UPON_LEAST_SQUARE_SOLUTION2D Sets up the equation system to solve for estimating a displacement field
%
% [a11, a12, a21, a22, b1, b2] = ...
%    estimate_displacement_based_upon_least_square_solution2d(solution, ...
%                                                             sizeImage, ...
%                                                             dk, ck, ...
%                                                             filterDirection, ...
%                                                             numberOfFilterDirections, ...
%                                                             t11, t12, t22)
%
% INPUT ARGUMENTS
% solution                      - Solution to use (0-4)
% sizeImage                     - Image size
% dk                            - Phase-difference for k different filter
%                                 directions
% ck                            - Corresponding certainties for k different
%                                 filter directions
% filterDirection               - k filter directions
% numberOfFilterDirections      - Number of filter directions
% t11                           - Tensor element (1,1)
% t12                           - Tensor element (1,2)
% t22                           - Tensor element (2,2)
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% a11 							- Components of local equation system A*d = b
% a12 							- Components of local equation system A*d = b
% a21 							- Components of local equation system A*d = b
% a22 							- Components of local equation system A*d = b
% b1 							- Components of local equation system A*d = b
% b2 							- Components of local equation system A*d = b

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
a21 = zeros(sizeImage);
a22 = zeros(sizeImage);
b1  = zeros(sizeImage);
b2  = zeros(sizeImage);

 switch solution
     case 0 %************************************************
         for k = 1:numberOfFilterDirections
             a11 = a11 + ck(:,:,k).*filterDirection{k}(1).^2;
             a12 = a12 + ck(:,:,k).*filterDirection{k}(1).*filterDirection{k}(2);
             a22 = a22 + ck(:,:,k).*filterDirection{k}(2).^2;

             b1  =  b1 + ck(:,:,k).*filterDirection{k}(1).*dk(:,:,k);
             b2  =  b2 + ck(:,:,k).*filterDirection{k}(2).*dk(:,:,k);
         end
         a21 = a12;
     case 1 %************************************************
         for k = 1:numberOfFilterDirections
             a11 = a11 + ck(:,:,k);
             a22 = a22 + ck(:,:,k);

             b1  =  b1 + ck(:,:,k).*dk(:,:,k).*filterDirection{k}(1);
             b2  =  b2 + ck(:,:,k).*dk(:,:,k).*filterDirection{k}(2);
         end
     case 2 %************************************************
         for k = 1:numberOfFilterDirections
             z1 = t11.*filterDirection{k}(1) + t12.*filterDirection{k}(2);
             z2 = t12.*filterDirection{k}(1) + t22.*filterDirection{k}(2);

             a11 = a11 + ck(:,:,k).*z1.^2;
             a12 = a12 + ck(:,:,k).*z1.*z2;
             a22 = a22 + ck(:,:,k).*z2.^2;

             b1  =  b1 + ck(:,:,k).*dk(:,:,k).*(filterDirection{k}(1).*z1.^2 + filterDirection{k}(2).*z1.*z2);
             b2  =  b2 + ck(:,:,k).*dk(:,:,k).*(filterDirection{k}(1).*z1.*z2 + filterDirection{k}(2).*z2.^2);
         end
         a21 = a12;
     case 3 %************************************************
         for k = 1:numberOfFilterDirections
             a11 = a11 + ck(:,:,k).*t11;
             a12 = a12 + ck(:,:,k).*t12;
             a22 = a22 + ck(:,:,k).*t22;

             b1  =  b1 + ck(:,:,k).*dk(:,:,k).*(filterDirection{k}(1)*t11 + (filterDirection{k}(2)*t12));
             b2  =  b2 + ck(:,:,k).*dk(:,:,k).*(filterDirection{k}(1)*t12 + (filterDirection{k}(2)*t22));
         end
         a21 = a12;
     case 4 %************************************************
         for k = 1:numberOfFilterDirections
             z1 = t11.*filterDirection{k}(1) + t12.*filterDirection{k}(2);
             z2 = t12.*filterDirection{k}(1) + t22.*filterDirection{k}(2);
             ztot = z1.^1 + z2.^2;

             a11 = a11 + ck(:,:,k).*filterDirection{k}(1).^2.*ztot;
             a12 = a12 + ck(:,:,k).*filterDirection{k}(1).*filterDirection{k}(2).*ztot;
             a22 = a22 + ck(:,:,k).*filterDirection{k}(2).^2.*ztot;

             b1  =  b1 + ck(:,:,k).*dk(:,:,k).*ztot.*filterDirection{k}(1);
             b2  =  b2 + ck(:,:,k).*dk(:,:,k).*ztot.*filterDirection{k}(2);
         end
         a21 = a12;
 end