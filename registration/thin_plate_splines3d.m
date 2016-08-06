function [output, displacementField] = thin_plate_splines3d(input, outputSize, ...
    inputLandmarks, outputLandmarks)
% THIN_PLATE_SPLINES3D Deforms input volume according to given landmarks and thin plate splines
%
% [output] = thin_plate_splines3d(input, outputSize, ...
%     inputLandmarks, outputLandmarks)
% 
% INPUT ARGUMENTS
% input             - Input volume, data to deformed
% outputSize        - Size of output volume
% inputLandmarks    - Landmarks in input volume
% outputLandmarks   - Landmarks in output volume
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% output            - Deformed data
% displacementField - Displacement field to deform input to output

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

%% Initialization
% Number of landmark points
NPs = size(inputLandmarks,1); 

% Size
height = outputSize(1);
width = outputSize(2);
depth = outputSize(3);

% Landmarks in input
Xi = inputLandmarks(:,1)';
Yi = inputLandmarks(:,2)';
Zi = inputLandmarks(:,3)';

% Landmarks in output (homologous)
Xo = outputLandmarks(:,1)';
Yo = outputLandmarks(:,2)';
Zo = outputLandmarks(:,3)';

%% Algebra of Thin-plate splines

% Compute thin-plate spline mapping [W|a1 ax ay az] using landmarks
[L]=computeWl3d(Xi, Yi, Zi, NPs);
% Y = ( V| 0 0 0)'   where V = [G] where G is landmark homologous (nx2) ; Y is col vector of length (n+3)
Y = [Xo(:) Yo(:) Zo(:); zeros(4,3)]; 
W = L\Y; % (W|a1 ax ay az)' = inv(L)*Y

% Thin-plate spline mapping (Map all points in the plane)
[Xw, Yw, Zw] = tpsMap(W, width, height, depth, Xi, Yi, Zi, NPs);

output = ba_interp3(input,Xw,Yw,Zw,'linear');
displacementField = {Xw,Yw,Zw};
end

%% [L] = [[K P];[P' 0]]
% np - number of landmark points
% (xp, yp, zp) - coordinate of landmark points
function [wL]=computeWl3d(xp, yp, zp, np)

rXp = repmat(xp(:),1,np); % 1xNp to NpxNp
rYp = repmat(yp(:),1,np); % 1xNp to NpxNp
rZp = repmat(zp(:),1,np); % 1xNp to NpxNp

wR = sqrt((rXp-rXp').^2 + (rYp-rYp').^2 + (rZp-rZp').^2); % compute r(i,j)

wK = radialBasis(wR); % compute [K] with elements U(r)=r^2 * log (r^2)
wP = [ones(np,1) xp(:) yp(:) zp(:)]; % [P] = [1 xp' yp' zp'] where (xp',yp',zp') are n landmark points (nx3)
wL = [wK wP;wP' zeros(4,4)]; % [L] = [[K P];[P' 0]]

end

%%
function [Xw, Yw, Zw] = tpsMap(W, width, height, depth, Xi, Yi, Zi, NPs)
    [xGrid,yGrid,zGrid] = meshgrid(1:width,1:height,1:depth); % HxW
    xGrid = xGrid(:)'; % convert to 1D array by reading columnwise (NWs=H*W*D)
    yGrid = yGrid(:)'; % convert to 1D array (NWs)
    zGrid = zGrid(:)'; % convert to 1D array (NWs)
    NWs = length(xGrid); % total number of points in the plane
    
    % All points in plane
    rX = repmat(xGrid,NPs,1); % Np x NWs
    rY = repmat(yGrid,NPs,1); % Np x NWs
    rZ = repmat(zGrid,NPs,1); % Np x NWs
    
    % Landmark points
    rxp = repmat(Xi(:),1,NWs); % 1xNp to Np x NWs
    ryp = repmat(Yi(:),1,NWs); % 1xNp to Np x NWs
    rzp = repmat(Zi(:),1,NWs); % 1xNp to Np x NWs
    
    % Mapping Algebra
    wR = sqrt((rxp-rX).^2 + (ryp-rY).^2 + (rzp-rZ).^2); % distance measure r(i,j)=|Pi-(x,y)|
    
    wK = radialBasis(wR); % compute [K] with elements U(r)=r^2 * log (r^2)
    wP = [ones(NWs,1) xGrid(:) yGrid(:) zGrid(:)]'; % [P] = [1 x' y'] where (x',y') are n landmark points (nx2)
    wL = [wK;wP]'; % [L] = [[K P];[P' 0]]
    
    Xw  = reshape(wL*W(:,1),[height width depth]); % [Pw] = [L]*[W]
    Yw  = reshape(wL*W(:,2),[height width depth]); % [Pw] = [L]*[W]
    Zw  = reshape(wL*W(:,3),[height width depth]); % [Pw] = [L]*[W]
end

%% k=(r^2) * log(r^2)
function [ko]=radialBasis(ri)

r1i = ri;
r1i(ri==0) = realmin; % Avoid log(0)=inf
ko = 2*(ri.^2).*log(r1i);

end