function [transformationMatrix] = ...
    compute_transformation_matrix_from_coordinates(newCoord,oldCoord,varargin)
% COMPUTE_TRANSFORMATION_MATRIX_FROM_COORDINATES
%
% Computes a transformation matrix that that transforms old coordinates to new
% coordinates
%
% [transformationMatrix] = ...
%     compute_transformation_matrix_from_coordinates(newCoord,oldCoord)
%
% INPUT ARGUMENTS
% newCoord
% oldCoord
%
% OPTIONAL INPUT ARGUMENTS
% 'dataSize'
%
% OUTPUT ARGUMENTS
% transformationMatrix

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

%% Default parameters
dataSize = [];

% Overwrites default parameters
for k=1:2:length(varargin)
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

%%
dims = size(newCoord,2);

if isempty(dataSize)
    % Do nothing
else
    % Compensate coordinates for data size, i.e. to make sure that the estimated
    % transformation has its center of rotation in the middle of the data volume
    dataSize(1:2) = fliplr(dataSize(1:2));
    newCoord = bsxfun(@minus,newCoord,dataSize/2);
    oldCoord = bsxfun(@minus,oldCoord,dataSize/2);
end

newCoord = [newCoord ones(size(newCoord,1),1)];
oldCoord = [oldCoord ones(size(oldCoord,1),1)];

T = pinv(oldCoord) * newCoord;

oldCoordTransformed = oldCoord * T;
if dims == 2
    transformationMatrix = transpose(T(:,1:2));
elseif dims == 3
    transformationMatrix = transpose(T(:,1:3));
    figure(23456)
    scatter3(oldCoord(:,1),oldCoord(:,2),oldCoord(:,3),'+r')
    hold on
    scatter3(oldCoordTransformed(:,1),oldCoordTransformed(:,2),oldCoordTransformed(:,3),'+g')
    scatter3(newCoord(:,1),newCoord(:,2),newCoord(:,3),'+b')
    hold off
    axis image
    xlabel('x-axis')
    ylabel('y-axis')
    zlabel('z-axis')
else
    error('Only 2D or 3D is currently supported')
end
