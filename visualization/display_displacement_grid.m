function display_displacement_grid(displacement,figNo,title)
% DISPLAY_DISPLACEMENT_GRID Displays the deformed grid
%
% display_displacement_grid(displacement,figNo,title)
%
% INPUT ARGUMENTS
% displacement  - Displacement field to display
% figNo     - Figure number to use
% title     - Figure title
%
% OPTIONAL INPUT ARGUMENTS
% N/A
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

figure(figNo)
imshow(zeros(size(displacement{1})),'InitialMagnification',100,'Border','tight')
figure_title(title)
hold on
plot_grid(displacement,'spacing',[10 10],'color','y','linewidth',2);
hold off
drawnow
