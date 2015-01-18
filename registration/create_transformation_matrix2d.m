function [transformationMatrix, rotationMatrix, shearMatrix, scaleMatrix] = ...
    create_transformation_matrix2d(varargin)
% CREATE_TRANSFORMATION_MATRIX2D Create a 2D transformation matrix
%
% [transformationMatrix rotationMatrix shearMatrix scaleMatrix] = ...
%   create_transformation_matrix2d()
% 
% INPUT ARGUMENTS
% N/A
%
% OPTIONAL INPUT ARGUMENTS
% 'phiZ'            - Rotation angle around z-axis
% 'shearXY'         - Shear XY
% 'shearYX'         - Shear YX
% 'scaleX'          - Scale along x-axis
% 'scaleY'          - Scale along y-axis
% 'translationX'    - Translation along x-axis
% 'translationY'    - Translation along y-axis
%
% OUTPUT ARGUMENTS
% transformationMatrix  - Total transformation matrix
% rotationMatrix        - Rotation matrix
% shearMatrix           - Shear matrix
% scaleMatrix           - Scale matrix

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

phiZ = 0;
shearXY = 0;
shearYX = 0;
scaleX = 1;           
scaleY = 1;
translationX = 0;
translationY = 0;

% Overwrites default parameter
for k=1:2:length(varargin)      
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

rotationMatrix = [cos(phiZ) -sin(phiZ); 
                  sin(phiZ) cos(phiZ)];

shearMatrix = [1       shearXY;
               shearYX 1      ];

scaleMatrix = [scaleX 0     ;
               0      scaleY];
           
transformationMatrix = zeros(2,3);
transformationMatrix(1:2,1:2) = rotationMatrix*shearMatrix*scaleMatrix;
transformationMatrix(:,3) = [translationX translationY]';
