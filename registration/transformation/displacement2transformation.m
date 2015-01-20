function [transformation] = displacement2transformation(displacement, imSize)
% DISPLACEMENT2TRANSFORMATION Converts a displacement field to a transformation matrix
%
% [transformation] = displacement2transformation(displacement, imSize)
%
% INPUT ARGUMENTS
% displacement          - Resulting displacement field
% imSize                - Image size
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% transformation  		- Transformation to convert

%
% Rotationcenter is in the center of the data

% Copyright (c) 2014 Daniel Forsberg
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

dims = length(imSize);

switch dims
    case 2
        sy = imSize(1);
        sx = imSize(2);
        [x,y] = meshgrid(linspace(-sx/2+1/2,sx/2-1/2,sx),linspace(-sy/2+1/2,sy/2-1/2,sy));
        
        xNew = vec(x + displacement{1});
        yNew = vec(y + displacement{2});
        
        newCoord = vec([xNew'; yNew']);
        
        x = x(:);
        y = y(:);
        oldCoord = [vec([x'; zeros(size(x'))]) vec([y'; zeros(size(x'))])...
            vec([ones(size(x')); zeros(size(x'))]) vec([zeros(size(x')) x'])...
            vec([zeros(size(x')); y']) vec([zeros(size(x')); ones(size(x'))])];
        
        T = pinv(oldCoord) * newCoord;
        
        transformation = zeros(2,3);
        
        transformation(1,1) = T(1);
        transformation(1,2) = T(2);
        transformation(1,3) = T(3);
        transformation(2,1) = T(4);
        transformation(2,2) = T(5);
        transformation(2,3) = T(6);
    case 3
        
    otherwise
        error('Only 2D or 3D is currently supported')
end
