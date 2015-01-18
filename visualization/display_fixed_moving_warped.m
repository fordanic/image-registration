function display_fixed_moving_warped(fixed,moving,warped,gamma)
% DISPLAY_FIXED_MOVING_WARPED Display three volumes in three planes
%
% display_fixed_moving_warped(fixed,moving,warped,gamma)
%
% INPUT ARGUMENTS
% fixed             - Fixed volume
% moving            - Moving volume
% warped            - Warped volume
% gamma             - Gamma value to apply
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

%**************** Display fixed ******************************
display_image(fixed,1,'fixed','gamma',gamma)

%**************** Display moving ******************************
display_image(moving,2,'moving','gamma',gamma)

%**************** Display warped ******************************
display_image(warped,3,'deformed','gamma',gamma)

%**************** Display (moving - fixed) *****************
display_images(moving,fixed,4,'moving-fixed','gamma',gamma)

%**************** Display (warped - fixed) *****************
display_images(warped,fixed,5,'deformed-fixed','gamma',gamma)
