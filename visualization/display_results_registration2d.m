function display_results_registration2d(moving,fixed,...
    displacementField,originalSize,gamma)
% DISPLAY_RESULTS_REGISTRATION2D Displays registration results for 2D
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
displacementField = resampler(displacementField, originalSize, 'relativeValues', true);
deformed = deformation(moving,displacementField,'linear');

%**************** Display fixed ******************************
display_image(fixed,1,'fixed','gamma',gamma);

%**************** Display moving ******************************
display_image(moving,2,'moving','gamma',gamma);

%**************** Display deformed ******************************
display_image(deformed,3,'deformed','gamma',gamma);

%**************** Display (moving - fixed) *****************
display_images(moving,fixed,4,'moving-fixed','gamma',gamma)

%**************** Display (deformed - fixed) *****************
display_images(deformed,fixed,5,'deformed-fixed','gamma',gamma)

%**************** Display displacement ******************************
display_displacement_grid(displacementField,6,'grid');

%************ Display accumulated motion field ********
if exist('display_motion_field','file')
    display_motion_field(displacementField,7,'Acc motion field',false)
end

%************ Display det of acc motion field **********
if exist('display_determinant_of_motion_field','file')
    display_determinant_of_motion_field(displacementField,8,'Det of jacobian of acc motion field',gamma)
end
