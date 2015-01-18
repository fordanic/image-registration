function [Y, X, incrementalCertainty] = ...
    estimate_displacement_based_upon_vector_summation2d(solution,...
    filterDirection,...
    sizeImage,...
    numberOfFilterDirections,...
    dk,ck,...
    t11,t22,t12)
% ESTIMATE_DISPLACEMENT_BASED_UPON_VECTOR_SUMMATION2D Estimates a displacement field based upon vector summation
%
% [Y, X, incrementalCertainty] = ...
%    estimate_displacement_based_upon_vector_summation2d(solution,...
%                                                        filterDirection,...
%                                                        sizeImage,...
%                                                        numberOfFilterDirections,...
%                                                        dk,ck,...
%                                                        t11,t22,t12)
%
% INPUT ARGUMENTS
% solution                      - Solution to use (6-7)
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
% Y	 							- Displacement along Y
% X 							- Displacement along X
% incrementalCertainty 			- Certainty of displacement

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

X = zeros(sizeImage);
Y = zeros(sizeImage);
incrementalCertainty = zeros(sizeImage);

switch solution
    case 6 %************************************************
        for k = 1:numberOfFilterDirections
            ckprim(:,:,k) = ck(:,:,k).*((filterDirection{k}(1).*t11 + filterDirection{k}(2).*t12).*filterDirection{k}(1)...
                +(filterDirection{k}(1).*t12 + filterDirection{k}(2).*t22).*filterDirection{k}(2));
        end
        
        for k = 1:numberOfFilterDirections
            X = X+ckprim(:,:,k).*dk(:,:,k)*filterDirection{k}(2);
            Y = Y+ckprim(:,:,k).*dk(:,:,k)*filterDirection{k}(1);
            incrementalCertainty = incrementalCertainty+ckprim(:,:,k);
        end
    case 7 %************************************************
        for k = 1:numberOfFilterDirections
            X = X+ck(:,:,k).*dk(:,:,k)*filterDirection{k}(2);
            Y = Y+ck(:,:,k).*dk(:,:,k)*filterDirection{k}(1);
            incrementalCertainty = incrementalCertainty+ck(:,:,k);
        end
end

incrementalCertainty = incrementalCertainty + eps;
X = X./incrementalCertainty;
Y = Y./incrementalCertainty;