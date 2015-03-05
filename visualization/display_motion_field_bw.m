function display_motion_field_bw(ima,motionField,figNo,title,removeMean)
% DISPLAY_MOTION_FIELD_BW Displays a motion field in BW
%
% display_motion_field_bw(ima,motionField,figNo,title,gamma,removeMean)
%
% INPUT ARGUMENTS
% motionField           - Motionfield to display
% figNo                 - Figure number to display motion field in
% title                 - Figure title
% gamma                 - gamma used in mygopimage
% removeMean            - Set to true to remove mean before displaying
%                         motion field
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

gamma = 1;
[sizeY sizeX] = size(motionField{1});
delta = 15;
xc = 1:delta:sizeX; 
yc = [1:delta:sizeY]';

%***************** Display motion field ****************
motionFieldPlot = motionField{1} + 1i*motionField{2};
motionFieldPlotRe = real(motionFieldPlot);
motionFieldPlotIm = imag(motionFieldPlot);
if removeMean
    motionFieldPlotRe = motionFieldPlotRe - mean(motionFieldPlotRe(:));
    motionFieldPlotIm = motionFieldPlotIm - mean(motionFieldPlotIm(:));
end
norm = sqrt(motionFieldPlotRe.^2 + motionFieldPlotIm.^2 + eps);
motionFieldPlotRe = motionFieldPlotRe./(norm.^gamma);
motionFieldPlotIm = motionFieldPlotIm./(norm.^gamma);
figure(figNo), gimage(ima/2,max(ima(:)),gamma,1), hold on
h = myquiver(xc, yc, motionFieldPlotRe(yc,xc), motionFieldPlotIm(yc,xc),1,'w');
set(gcf, 'name', title)
set(h,'linewidth',1)
hold off
drawnow
