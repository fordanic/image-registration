function display_displacement_grid(displacement,figNo,title,varargin)
% DISPLAY_DISPLACEMENT_GRID Displays the deformed grid
%
% display_displacement_grid(displacement,figNo,title)
%
% INPUT ARGUMENTS
% displacement  - Displacement field to display
% figNo         - Figure number to use
% title         - Figure title
%
% OPTIONAL INPUT ARGUMENTS
% 'spacing'     - spacing of grid
% 'color'       - Color of grid
% 'linewidth'   - Line width of grid
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

%% Default parameters
% Setup default parameter
spacing = [10 10];
color = 'b';
linewidth = 2;

% Overwrites default parameter
for k=1:2:length(varargin)
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

%%
figure(figNo)
imshow(zeros(size(displacement{1})),'InitialMagnification',100,'Border','tight')
figure_title(title)
hold on
plot_grid(displacement,'spacing',spacing,'color',color,'linewidth',linewidth);
hold off
drawnow
