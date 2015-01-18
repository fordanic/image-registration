function setup_image_registration_repository()
% SETUP_IMAGE_REGISTRATION_REPOSITORY Call this function to setup and compile necessary files
%
% setup_image_registration_repository

% Copyright (c) 2014 Daniel Forsberg
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

%% Compile mex files
[repositoryPath,~,~] = fileparts(mfilename('fullpath'));

compile_mex_files([repositoryPath,filesep,'external',filesep,'misc']);

compile_mex_files([repositoryPath,filesep,'external',filesep,'spatial_domain_toolbox'],...
    'filesToExclude',{'antigradient.c'});

compile_mex_files([repositoryPath,filesep,'registration',filesep,'polynomial-expansion'],...
    'fileExtension','cpp');

%% Create various filters

% Quadrature filters for linear phase-difference registration
if ~exist('quadratureFiltersLinearRegistration2D.mat')
    create_quadrature_filters_for_phase_difference_registration2d();
end
if ~exist('quadratureFiltersLinearRegistration3D.mat')
    create_quadrature_filters_for_phase_difference_registration3d();
end

% quadrature filter for morphon registration
if ~exist('quadratureFiltersForMorphonRegistration2D.mat')
    create_quadrature_filters_for_morphon_registration2d();
end
if ~exist('quadratureFiltersForMorphonRegistration3D.mat')
    create_quadrature_filters_for_morphon_registration3d();
end