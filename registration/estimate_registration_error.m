function RE = estimate_registration_error(fixedCoordinates,movingCoordinates,displacementField,spacing)
% ESTIMATE_REGISTRATION_ERROR Estimates the target registration error
% 
% RE = estimate_registration_error(...
%   fixedCoordinates,movingCoordinates,displacementField,spacing)
%
% INPUT ARGUMENTS
% fixedCoordinates      - Coordinates in fixed image
% movingCoordinates     - Coordinates in moving image
% displacementField     - Displacement field, pointing from fixed to moving
% spacing               - Pixel spacing
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% RE                    - Registration error

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

dims = length(displacementField);
deformedFixedCoordinates = fixedCoordinates;
switch dims
    case 1
        xDisp = interp3(displacementField{1},fixedCoordinates(:,1));
        deformedFixedCoordinates(:,1) = fixedCoordinates(:,1) + xDisp;
        RE = deformedFixedCoordinates - movingCoordinates;
        RE = sqrt((RE(:,1)*spacing(1)).^2);
    case 2
        xDisp = interp3(displacementField{1},fixedCoordinates(:,1),fixedCoordinates(:,2));
        yDisp = interp3(displacementField{2},fixedCoordinates(:,1),fixedCoordinates(:,2));
        deformedFixedCoordinates(:,1) = fixedCoordinates(:,1) + xDisp;
        deformedFixedCoordinates(:,2) = fixedCoordinates(:,2) + yDisp;
        RE = deformedFixedCoordinates - movingCoordinates;
        RE = sqrt((RE(:,1).^spacing(1)).^2 + (RE(:,2)*spacing(2)).^2);
    case 3
        xDisp = interp3(displacementField{1},fixedCoordinates(:,1),fixedCoordinates(:,2),fixedCoordinates(:,3));
        yDisp = interp3(displacementField{2},fixedCoordinates(:,1),fixedCoordinates(:,2),fixedCoordinates(:,3));
        zDisp = interp3(displacementField{3},fixedCoordinates(:,1),fixedCoordinates(:,2),fixedCoordinates(:,3));
        deformedFixedCoordinates(:,1) = fixedCoordinates(:,1) + xDisp;
        deformedFixedCoordinates(:,2) = fixedCoordinates(:,2) + yDisp;
        deformedFixedCoordinates(:,3) = fixedCoordinates(:,3) + zDisp;
        RE = deformedFixedCoordinates - movingCoordinates;
        RE = sqrt((RE(:,1)*spacing(1)).^2 + (RE(:,2)*spacing(2)).^2 + (RE(:,3)*spacing(3)).^2);
    otherwise
        error('Not implemented for dimensions higher than three.')
end