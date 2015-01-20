function transformationMatrix = create_transformation_matrix3d(varargin)
% CREATE_TRANSFORMATION_MATRIX3D Create a 3D transformation matrix
%
% [transformationMatrix rotationMatrix shearMatrix scaleMatrix] = ...
%   create_transformation_matrix3d()
% 
% INPUT ARGUMENTS
% N/A
%
% Optional input arguments
% 'phiX'            - Rotation angle around x-axis
% 'phiY'            - Rotation angle around y-axis
% 'phiZ'            - Rotation angle around z-axis
% 'shearXY'         - Shear XY
% 'shearXZ'         - Shear XZ
% 'shearYX'         - Shear YX
% 'shearYZ'         - Shear YZ
% 'shearZX'         - Shear ZX
% 'shearZY'         - Shear ZY
% 'scaleX'          - Scale along x-axis
% 'scaleY'          - Scale along y-axis
% 'scaleZ'          - Scale along z-axis
% 'translationX'    - Translation along x-axis
% 'translationY'    - Translation along y-axis
% 'translationZ'    - Translation along z-axis
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

phiX = 0;
phiY = 0;
phiZ = 0;
shearXY = 0;
shearXZ = 0;
shearYX = 0;
shearYZ = 0;
shearZX = 0;
shearZY = 0;
scaleX = 1;           
scaleY = 1;
scaleZ = 1;
translationX = 0;
translationY = 0;
translationZ = 0;

for k=1:2:length(varargin),         % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

rotationX = [1 0         0; 
             0 cos(phiX) -sin(phiX); 
             0 sin(phiX) cos(phiX)];
rotationY = [cos(phiY)  0 sin(phiY); 
             0          1 0; 
             -sin(phiY) 0 cos(phiY)];
rotationZ = [cos(phiZ) -sin(phiZ) 0; 
             sin(phiZ) cos(phiZ)  0; 
             0         0          1];

shearMatrix = [1       shearXY shearXZ;
               shearYX 1       shearYZ;
               shearZX shearZY 1];

scaleMatrix = [scaleX 0      0;
               0      scaleY 0;
               0      0      scaleZ];
           

transformationMatrix = zeros(3,4);
transformationMatrix(1:3,1:3) = rotationX*rotationY*rotationZ*shearMatrix*scaleMatrix;
transformationMatrix(:,4) = [translationX translationY translationZ]';
