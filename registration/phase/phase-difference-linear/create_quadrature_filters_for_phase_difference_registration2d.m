function create_quadrature_filters_for_phase_difference_registration2d(varargin)
% CREATE_QUADRATURE_FILTERS_FOR_PHASE_DIFFERENCE_REGISTRATION2D Creates 2D quadrature filters for usage in linear phase-difference registration
%
% create_quadrature_filters_for_phase_difference_registration2d()
%
% INPUT ARGUMENTS
% N/A
%
% OPTIONAL INPUT ARGUMENTS
% See "Advanced filter design" by Knutsson et al for detailed explanation
% of the parameters used. Here, only default values are listed.
%
%                       Default values
% 'filterSize'          - 11
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

filePath = fileparts(mfilename('fullpath'));
currentPath = pwd;

cd(filePath);

%% Set deafult parameters

% Filter parameters
filterSize = 11;
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

fprintf('Optimizing quadrature filters for linear phase difference registration in 2d\n');

% Sizes
spatialSize = [filterSize filterSize];
frequencySize = 2*spatialSize+1;

% Frequency weights
Fw = F_weightgensnr(frequencySize,frequency_rexp,cosexp,SNRest,DCamp);

% Spatial weights
fw0 = goodsw(spatialSize,spatial_rexp);
%wamesh(fw0)
fw_amp = 30;
fw = fw_amp.*fw0;

% Spatial ideal filter
fi = wa(spatialSize,0);
fi = putorigo(fi,1);

% Spatial mask
fm = wa(spatialSize,0);
fm = putdata(fi,ones(spatialSize));

% Filter directions
dir{1} = [0 1]'; % x
dir{2} = [1 0]'; % y

% Frequency ideal filters
for k = 1 : 2
    Fi{k} = quadrature(frequencySize,u0,B,dir{k});
end

% Optimize the quadrature filters
for k = 1 : 2
    [f{k},F{k}] = krnopt(Fi{k},Fw,fm,fi,fw);
end

filterr(F{1},Fi{1},Fw,fi,fw,1);

f1 = getdata(f{1});
f2 = getdata(f{2});

save quadratureFiltersLinearRegistration2D f1 f2

cd(currentPath)