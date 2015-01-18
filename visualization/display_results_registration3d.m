function display_results_registration3d(moving,fixed,...
    displacementField,originalSize,gamma)
% DISPLAY_RESULTS_REGISTRATION3D Displays registration results for 3D
%
% display_results_registration3d(moving,fixed,displacementField,originalSize,gamma)
%
% INPUT ARGUMENTS
% moving                - Moving image
% fixed                 - Fixed image
% displacementField     - Displacement field
% originalSize          - Original data size, data is resampled to this
%                         size before display
% gamma                 - Gamma to use for display
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

moving = moving/max(moving(:));
fixed = fixed/max(fixed(:));
displacementField = resampler(displacementField, originalSize, 'relativeValues', 1);
centerPosition = floor(size(moving) / 2);

deformed = deformation(moving,displacementField,'linear');
deformed = deformed/max(deformed(:));

%**************** Display fixed, moving and warped ******************************
display_fixed_moving_warped(fixed,moving,deformed,gamma);

%**************** Display displacement field ******************************
disp{1} = transpose(squeeze(displacementField{1}(:,:,centerPosition(3))));
disp{2} = transpose(squeeze(displacementField{2}(:,:,centerPosition(3))));
display_displacement_grid(disp,10,'grid YX');

disp{1} = transpose(squeeze(displacementField{2}(:,centerPosition(2),:)));
disp{2} = transpose(squeeze(displacementField{3}(:,centerPosition(2),:)));
display_displacement_grid(disp,11,'grid ZY');

disp{1} = transpose(squeeze(displacementField{1}(centerPosition(1),:,:)));
disp{2} = transpose(squeeze(displacementField{3}(centerPosition(1),:,:)));
display_displacement_grid(disp,12,'grid ZX');
