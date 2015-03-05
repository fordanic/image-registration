function display_motion_field(motionField,figNo,title,removeMean,varargin)
% DISPLAY_MOTION_FIELD Displays a motion field
%
% display_motion_field(motionField,figNo,title,removeMean)
%
% INPUT ARGUMENTS
% motionField           - Motionfield to display
% figNo                 - Figure number to display motion field in
% title                 - Figure title
% removeMean            - Set to true to remove mean before displaying
%                         motion field
%
% OPTIONAL INPUT ARGUMENTS
% 'scaling'             - Scaling of quiver
% 'color'               - Color of quiver
% 'delta'               - Spacing for quiver
%
% OUTPUT ARGUMENTS
% N/A

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

% Setup default parameter
scaling = 1;
color = 'b';
delta = 15;

% Overwrites default parameter
for k=1:2:length(varargin)
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

[sizeY sizeX] = size(motionField{1});
xc = 1:delta:sizeX; 
yc = [1:delta:sizeY]';

%***************** Display motion field ****************
motionField{1} = motionField{1};
motionField{2} = motionField{2};
if removeMean
    motionField{1} = motionField{1} - mean(motionField{1}(:));
    motionField{2} = motionField{2} - mean(motionField{2}(:));
end
display_image(ones(size(motionField{1})),figNo,title);
hold on
h = quiver(xc, yc, motionField{1}(yc,xc), motionField{2}(yc,xc),scaling,color);
set(h,'linewidth',1.0)
hold off
drawnow
