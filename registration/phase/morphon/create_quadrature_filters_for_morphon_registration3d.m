function create_quadrature_filters_for_morphon_registration3d(varargin)
% CREATE_QUADRATURE_FILTERS_FOR_MORPHON_REGISTRATION3D Creates 3D quadrature filters for usage in Morphon registration
%
% create_quadrature_filters_for_morphon_registration3d()
%
% INPUT ARGUMENTS
% N/A
%
% OPTIONAL INPUT ARGUMENTS
% See "Advanced filter design" by Knutsson et al for detailed explanation
% of the parameters used. Here, only default values are listed.
%
%                       Default values
% 'filterSize'          - 9
% 'spatial_rexp'        - 2
% 'frequency_rexp'      - -1
% 'cosexp'              - 1
% 'SNRest'              - 30
% 'DCamp'               - 10000
% 'u0'                  - pi/3
% 'B'                   - 2.0
%
% OUTPUT ARGUMENTS
% N/A
%
% See also F_weightgensnr, goodsw, quadrature, krnopt

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

filePath = which('create_quadrature_filters_for_morphon_registration3d.m');
filePath = filePath(1:strfind(filePath,'create_quadrature_filters_for_morphon_registration3d')-1);
currentPath = pwd;

cd(filePath);

%% Set defaults parameters

% Filter parameters
filterSize = 9;
spatial_rexp = 2;
frequency_rexp = -1;
cosexp = 1;
SNRest = 30;
DCamp = 10000;

if ~exist('F_weightgensnr.m','file') || ~exist('goodsw.m','file') || ...
        ~exist('quadrature.m','file') || ~exist('krnopt.m','file')
    error(['You appear to be missing vital functions for optimizing filters. ',...
        'Please download: https://www.imt.liu.se/edu/courses/TBMI02/code/kerngen.zip'])
end

% Center frequency and bandwidth (in octaves)
u0 = pi/3;
B = 2.0;

%% Overwrites default parameter
for k=1:2:length(varargin)
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

fprintf('Optimizing quadrature filters for Morphon registration in 3d\n');

%*********************************************************
% Sizes
%*********************************************************
spatialSize = [filterSize filterSize filterSize];
frequencySize = 2*spatialSize+1;

%*****************************************************
% Spatial mask, fm
%*****************************************************
fm = ones(wa(spatialSize,0));          

%*****************************************************
% Spatial ideal filter, fi.
%*****************************************************
fi = wa(spatialSize, 0);
fi = putorigo(fi, 1);

%*****************************************************
% Spatial locality  weight, fw0, sum(fw0.^2) = 1
%*****************************************************
fw0 = goodsw(spatialSize, spatial_rexp);

%*****************************************************
% Frequency ideal filter, Fi
% same radial function as grad but cos^2
% as angular function
%*****************************************************

a = 2;
b = 1 + sqrt(5);
c = 1/sqrt(10 + 2*sqrt(5));
filterDirection{1} = c * [ 0  a  b];
filterDirection{2} = c * [ 0 -a  b];
filterDirection{3} = c * [ a  b  0];
filterDirection{4} = c * [-a  b  0];
filterDirection{5} = c * [ b  0  a];
filterDirection{6} = c * [ b  0 -a];

for k = 1 : length(filterDirection)
    Fi{k} = quadrature(frequencySize,u0,B,filterDirection{k});
end

    
%*****************************************************
% Frequency weight function, Fw.
% sum(Fw.^2) = 1 (excluding the origin)
%*****************************************************
Fw = F_weightgensnr(frequencySize, frequency_rexp, cosexp, SNRest, DCamp);

%*****************************************************
% optimize
%*****************************************************
fw_amp = 30;

fw = fw_amp.*fw0;
for k = 1 : length(filterDirection)
    fprintf('Optimizing quadrature filter %d\r',k)
    [f{k}, F{k}] = krnopt(Fi{k}, Fw, fm, fi, fw);
end

fprintf('                              \r',k)


%*****************************************************
% Show result
%*****************************************************
fprintf('                              Quad filter{1}\n')
filterr(F{1}, Fi{1}, Fw, f{1}, fw, 1);

fprintf('                              Quad filter{2}\n')
filterr(F{2}, Fi{2}, Fw, f{2}, fw, 1);

%*****************************************************
% Create direction tensor M associated with filters
%*****************************************************
for k = 1 : length(filterDirection)
    M = 5/4.*filterDirection{k}'*filterDirection{k} - 1/4.*eye(3);
    m11{k} = M(1,1);
    m12{k} = M(1,2);
    m13{k} = M(1,3);
    m22{k} = M(2,2);
    m23{k} = M(2,3);
    m33{k} = M(3,3);
end

%*****************************************************
% Save created filters
%*****************************************************
for k = 1 : length(filterDirection)
    f{k} = getdata(f{k});
    F{k} = getdata(F{k});
end

save('quadratureFiltersForMorphonRegistration3D','f', 'filterDirection', 'm11', 'm12', 'm13', 'm22', 'm23', 'm33')

cd(currentPath)