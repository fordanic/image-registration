function [displacement] = transformation2displacement(transformationMatrix, imSize)
% TRANSFORMATION2DISPLACEMENT Converts a transformation matrix to a displacement field
%
% [displacement] = transformation2displacement(transformationMatrix, imSize)
%
% INPUT ARGUMENTS
% transformationMatrix  - Transformation to convert
% imSize                - Image size
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% displacement          - Resulting displacement field
%
% Rotationcenter is in the center of the data

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

dims = length(imSize);

switch dims
    case 2
        sy = imSize(1);
        sx = imSize(2);
        [x,y] = meshgrid(linspace(-sx/2+1/2,sx/2-1/2,sx),linspace(-sy/2+1/2,sy/2-1/2,sy));
        
        M = zeros([size(x),3]);
        M(:,:,1) = x;
        M(:,:,2) = y;
        M(:,:,3) = 1;
        M = permute(M,[3,1,2]);
        
        xi = reshape(transformationMatrix(1,:)*M(:,1:sx*sy),[sy sx]);
        yi = reshape(transformationMatrix(2,:)*M(:,1:sx*sy),[sy sx]);
        
        displacement = cell(2,1);
        displacement{1} = xi - x;
        displacement{2} = yi - y;
    case 3
        sy = imSize(1);
        sx = imSize(2);
        sz = imSize(3);
        
        [x,y,z] = meshgrid(linspace(-sx/2+1/2,sx/2-1/2,sx),...
            linspace(-sy/2+1/2,sy/2-1/2,sy),...
            linspace(-sz/2+1/2,sz/2-1/2,sz));
        M = zeros([size(x),4]);
        M(:,:,:,1) = x;
        M(:,:,:,2) = y;
        M(:,:,:,3) = z;
        M(:,:,:,4) = 1;
        M = permute(M,[4,1,2,3]);
        
        xi = reshape(transformationMatrix(1,:)*M(:,1:sx*sy*sz),[sy sx sz]);
        yi = reshape(transformationMatrix(2,:)*M(:,1:sx*sy*sz),[sy sx sz]);
        zi = reshape(transformationMatrix(3,:)*M(:,1:sx*sy*sz),[sy sx sz]);
        
        displacement = cell(3,1);
        displacement{1} = xi - x;
        displacement{2} = yi - y;
        displacement{3} = zi - z;     
    otherwise
        error('Only 2D or 3D is currently supported')
end
