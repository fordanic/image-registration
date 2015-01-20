function [displacementField] = create_nonrigid_displacement_field3d(fieldSize,varargin)
% CREATE_NON_RIGID_DISPLACEMENT_FIELD3D Creates a 3D non-rigid displacement field
%
% [displacementField] = create_nonrigid_displacement_field3d(fieldSize)
%
% INPUT ARGUMENTS
% fieldSize             - [rowSize colSize rodSize], assumed to be of even size
%
% OPTIONAL INPUT ARGUMENTS
% 'maxDisplacement'     - Maximum displacement of the created displacement
%                         field (10 pixels, default value)
%
% OUTPUT ARGUMENTS
% displacementField     - Created displacement field

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
maxDisplacement = 10;

% Overwrites default parameter
for k=1:2:length(varargin)
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

%%
displacementField = cell(3,1);

for k = 1 : 3
	displacementField{k} = rand(fieldSize);
	displacementFieldff = fftn(displacementField{k});
	displacementFieldff = fftshift(displacementFieldff);
	sz = size(displacementField{k});
	[x,y,z]=meshgrid(-sz(2)/2:sz(2)/2-1,-sz(1)/2:sz(1)/2-1,-sz(3)/2:sz(3)/2-1);
	w = 1./(sqrt(x.^2 + y.^2 + z.^2));
	w(sz(2)/2+1,sz(1)/2+1,sz(3)/2+1) = 1;
	displacementFieldffw = displacementFieldff.*w;
	displacementFieldffw = ifftshift(displacementFieldffw);
	displacementFieldlp = ifftn(displacementFieldffw);
	displacementField{k} = 20*(displacementFieldlp - mean(displacementFieldlp(:)));
	displacementField{k} = 64*averaging3d(displacementField{k}, 0, 15, 1.5);
end

% Estimate current max displacement
displacement = max(sqrt(displacementField{1}(:).^2 + ...
    displacementField{2}(:).^2 + displacementField{3}(:).^2));

% If larger than set max displacement, limit.
if displacement > maxDisplacement
    for k = 1 : 3
        displacementField{k} = displacementField{k}/(displacement/maxDisplacement);
    end
end