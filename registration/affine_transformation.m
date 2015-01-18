function [outImage, displacement] = affine_transformation(inImage, transformMatrix, varargin)
% AFFINE_TRANSFORMATION Applies an affine transform to the provided image
%
% [outImage, displacement] = affine_transformation(inImage, transformMatrix)
%
% Rotationcenter is in the center of the data
%
% INPUT ARGUMENTS
% inImage                   - Image to transform
% transformMatrix           - Transformation to apply
%
% OPTIONAL INPUT ARGUMENTS
% 'interpolation'           - Interpolation method
%                           - 'nearest', 'linear' (default), 'cubic'
% 'spacing'                 - Set coordinate spacing, default is [1 1 1]
% 'useImCenterAsRotCenter'  - Set to true to use center of image as 
%                             rotation center, default is true
%
% OUTPUT ARGUMENTS
% outImage                  - Transformed image
% displacement              - Resulting displacement field

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

if isempty(transformMatrix)
    outImage = inImage;
    return;
end

%% Default parameters

% Set default parameters

interpolation = 'linear';
spacing = [1 1 1];
useImCenterAsRotCenter = true;

% Overwrites default parameters
for k=1:2:length(varargin)
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

%%
dims = ndims(inImage);

switch dims
    case 2
        [sy, sx] = size(inImage);
        if useImCenterAsRotCenter
            [x,y] = meshgrid(...
                linspace(-sx/2+1/2,sx/2-1/2,sx),linspace(-sy/2+1/2,sy/2-1/2,sy));
        else
            [x,y] = meshgrid(...
                linspace(1,sx,sx),linspace(1,sy,sy));
        end
        x = x * spacing(1);
        y = y * spacing(2);
        
        M = zeros([size(x),3]);
        M(:,:,1) = x;
        M(:,:,2) = y;
        M(:,:,3) = 1;
        M = permute(M,[3,1,2]);
        
        xi = reshape(transformMatrix(1,:)*M(:,1:sx*sy),[sy sx]);
        yi = reshape(transformMatrix(2,:)*M(:,1:sx*sy),[sy sx]);
        outImage = interp2(x,y,inImage,xi,yi,interpolation);
        outImage(isnan(outImage)) = 0;
        
        if nargout == 2
            displacement = cell(2,1);
            displacement{1} = xi - x;
            displacement{2} = yi - y;
        end
        
    case 3
        [sy, sx, sz] = size(inImage);
        
        if useImCenterAsRotCenter
            [x,y,z] = meshgrid(linspace(-sx/2+1/2,sx/2-1/2,sx),...
                linspace(-sy/2+1/2,sy/2-1/2,sy),...
                linspace(-sz/2+1/2,sz/2-1/2,sz));
        else
            [x,y,z] = meshgrid(...
                linspace(1,sx,sx),linspace(1,sy,sy),linspace(1,sz,sz));
        end
        x = x * spacing(1);
        y = y * spacing(2);
        z = z * spacing(3);
        M = zeros([size(x),4]);
        M(:,:,:,1) = x;
        M(:,:,:,2) = y;
        M(:,:,:,3) = z;
        M(:,:,:,4) = 1;
        M = permute(M,[4,1,2,3]);
        
        xi = reshape(transformMatrix(1,:)*M(:,1:sx*sy*sz),[sy sx sz]);
        yi = reshape(transformMatrix(2,:)*M(:,1:sx*sy*sz),[sy sx sz]);
        zi = reshape(transformMatrix(3,:)*M(:,1:sx*sy*sz),[sy sx sz]);
        outImage = ba_interp3(x,y,z,inImage,xi,yi,zi,interpolation);
        
        if nargout == 2
            displacement = cell(3,1);
            displacement{1} = xi - x;
            displacement{2} = yi - y;
            displacement{3} = zi - z;
        end
        
    otherwise
        error('Only 2D/3D is currently supported for affine')
end
