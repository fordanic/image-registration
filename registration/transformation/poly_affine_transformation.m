function [deformed, displacement] = poly_affine_transformation(...
    im, transformationMatrices,weightFunctions,varargin)
% POLY_AFFINE_TRANSFORMATION Fuses multiple affine transformations using a log-euclidean framework
%
% [deformed displacement] = poly_affine_transformation(...
%   im,transformationMatrices,weightFunctions)
%
% INPUT ARGUMENTS
% im                        - Image to deform	
% transformationMatrices    - K transformation matrices to fuse
% weightFunctions           - Weight-functions associated with each
%                             transformation matrix
%
% OPTIONAL INPUT ARGUMENTS
% 'N'                       - Number of "Scaling and Squaring" to apply (6)
%
% OUTPUT ARGUMENTS
% deformed                  - Deformed image
% displacement              - Displacement field associated with the
%                             poly-affine transformation

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

% Set up default parameters
N = 6;
GPU = false;

% Overwrites default parameter
for k=1:2:length(varargin)
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

sz = size(im);
dims = ndims(im);

switch dims
    case 2
        %% Compute initial displacement
        [x, y] = meshgrid(linspace(-sz(2)/2+1/2,sz(2)/2-1/2,sz(2)),...
            linspace(-sz(1)/2+1/2,sz(1)/2-1/2,sz(1)));
        displacement = cell(dims,1);
        displacement{1} = zeros(size(x));
        displacement{2} = zeros(size(y));
        for k = 1 : length(transformationMatrices)
            T = expm(1/2^N * logm(transformationMatrices{k}));
            xi = T(1,1).*x + T(1,2).*y + T(1,3);
            yi = T(2,1).*x + T(2,2).*y + T(2,3);
            
            displacement{1} = displacement{1} + weightFunctions{k}.*(xi - x);
            displacement{2} = displacement{2} + weightFunctions{k}.*(yi - y);
        end
    case 3
        %% Compute initial displacement
        [x, y, z] = meshgrid(linspace(-sz(2)/2+1/2,sz(2)/2-1/2,sz(2)),...
            linspace(-sz(1)/2+1/2,sz(1)/2-1/2,sz(1)),...
            linspace(-sz(3)/2+1/2,sz(3)/2-1/2,sz(3)));
        displacement = cell(dims,1);
        displacement{1} = zeros(size(x));
        displacement{2} = zeros(size(y));
        displacement{3} = zeros(size(3));
        for k = 1 : length(transformationMatrices)
            T = expm(1/2^N * logm(transformationMatrices{k}));
            xi = T(1,1).*x + T(1,2).*y + T(1,3).*z + T(1,4);
            yi = T(2,1).*x + T(2,2).*y + T(2,3).*z + T(2,4);
            zi = T(3,1).*x + T(3,2).*y + T(3,3).*z + T(3,4);
            
            displacement{1} = displacement{1} + weightFunctions{k}.*(xi - x);
            displacement{2} = displacement{2} + weightFunctions{k}.*(yi - y);
            displacement{3} = displacement{3} + weightFunctions{k}.*(zi - z);
        end
    otherwise
        error('Only supported for 2D/3D data');
end

%% Perform N field compositions
if GPU
    [displacement{1}, displacement{2}, displacement{3}] = ...
        CUDA_field_composition3d(displacement{1}, displacement{2}, displacement{3}, N);
else
    newDisplacement = cell(dims,1);
    for l = 1 : N
        for k = 1 : dims
            newDisplacement{k} = deformation(displacement{k},displacement,'linear');
            newDisplacement{k} = newDisplacement{k}+displacement{k};
        end
        displacement = newDisplacement;
    end
end

deformed = deformation(im,displacement,'linear');
