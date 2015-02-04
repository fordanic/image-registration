# image-registration

A MATLAB library/toolbox providing access to image registration, 
initially developed for use with medical images.

# Copyright

Copyright (c) 2015 Daniel Forsberg
danne.forsberg@outlook.com

# About

This toolbox was written during my PhD at Linköping University, Sweden. 
However, although developed during that time, it doesn't contain so much
of the specifics that I researched during my PhD. Rather it contains a 
general framework for image registration in MATLAB. The toolbox supports 
the following features:
* 2D and 3D registration
* Symmetric registration
* Translation, rigid and affine registration (parametric)
* Non-rigid registration (non-parametric) with support for
** Fluid and elastic regularization (not really but often called so)
** Additive, compositive and diffeomorphic accumulation of the update field
** Registration in the Log-domain
* All this using either optical-flow based registration (demons in the 
non-rigid case), phase-difference based registration (the Morphon in the
non-rigid case) and polynomial expansion based registration.
* Use of certainty fields for controlling areas of influence during the 
registration process (only in the case of phase-difference based registration)

Specifics of my thesis are rather found in separate toolboxes, 
to be made available at my GitHub account.

# Setup

To use the code available in this repository, add the following 
lines to your startup.m file.

addpath('path-to-this-repository')

setup_image_registration_repository();

The last line will make sure to compile included mex-files and 
keep them up to date.

If you don't have a startup.m file, create it from startupsav.m.
Run the 'userpath' function to determine where MATLAB looks for 
the startup.m file.

Note that this library is dependent on my matlab-utilities repository,
available from https://github.com/fordanic/matlab-utilities and assumes
it to already be installed and available on the MATLAB path. If not 
available then you will be offered during setup to download needed
external dependencies.

Best way to test the setup of the toolbox is to run registration_example2d 
or registration_example3d. Type help registration_example2d
or help registration_example3d to see some of the parameters
that can be set. Use these examples as a starting point
when setting up to run your own registration.

# Content

This section provides a brief overview of the content in the 
different folders

external - This folder contains code from external sources

registration - This folder contains functions for image registration

test-data - This folder contains data for testing the registration.

visualization - This folder contains functions related to 
visualization of the registration process

# Adding folders

When adding a new folder in the base folder make sure to update 
this file.

# Coding standard

Basic guidelines for a simple coding standard are given in the document 
"Matlab Programming Style Guidelines.pdf", available in the
matlab-utilities repository.

# Adding m-files

To create a suitable file header for new M-files, use the function 
create_new_m_function, available in the matlab-utilities repository.

# Licensing

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Source code provided in this repository is generally released under 
the GNU GENERAL PUBLIC LICENSE Version 3, if not otherwise stated.
