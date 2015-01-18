function [v, regions, params, startpoints] = velocity_segment(T, options)
    
% [V, REGIONS, PARAMS, STARTPOINTS] = VELOCITY_SEGMENT(T, OPTIONS)
% 
% Compute simultaneous segmentation and velocity estimation from the
% spatiotemporal orientation tensor field T. The algorithm is described
% in section 6.4 of Gunnar Farnebäck's thesis "Polynomial Expansion for
% Orientation and Motion Estimation".
% 
% T
%   Orientation tensor field for one frame, an array of size
%   height x width x 3 x 3.
%
% OPTIONS [optional]
%   Struct array that may contain various parameters that affect the
%   algorithm. These are explained below. The default value is no options.
%
% V
%   Estimated velocity field of size height x width x 2.
% 
% REGIONS [optional]
%   Integer valued matrix of size height x width. Pixels with the same
%   number belong to the same region. Zero-valued pixels are unassigned.
%
% PARAMS [optional]
%   Cell array containing motion model parameters for each region.
%
% STARTPOINTS [optional]
%   2xM matrix, where M is the number of regions. Each column gives the
%   coordinate for one point in the corresponding region. This was used
%   as starting point in the region growing process.
%
%
% The following fields, each being a scalar value, may be specified in
% the OPTIONS parameter:
%
% OPTIONS.lambda -
%                 Relative cost of adding a new region compared to
%                 extending an existing region. A high lambda value
%                 leads to a few large regions while a small value
%                 gives a highly fragmented result. If good velocity
%                 estimates are most important a small value is usually
%                 best. See section 6.4.3 for more information.
%                 Default value is 0.1 (fairly small).
%
% OPTIONS.maxsize -
%                 Size a candidate region is (re)grown to and has
%                 when it is elevated to a real region. See section
%                 6.4.2 for more information, where it is called m_0.
%                 Default value is 500.
%                 
% OPTIONS.minsize -
%                 If a candidate region is so restricted by its
%                 neighbors that it cannot be (re)grown to this size,
%                 it is dicarded. This parameter is not discussed in the
%                 thesis, where it is assumed to be the same as maxsize.
%                 Default value is 500. If it turns out higher than maxsize,
%                 it is reset to maxsize.
%
% OPTIONS.kerneldist -
%                 Spacing of the starting points for candidate regions.
%                 A smaller value is slower while a larger value
%                 potentially leads to a worse segmentation. Default value
%                 is 4.
%                 
% OPTIONS.coverage -
%                 How large part of the image to segment. Default value 
%                 is 1.0 (i.e. the entire image).
%                 
% OPTIONS.parameter_recomputation -
%                 Whether the motion model parameters are recomputed
%                 at the end of the algorithm. Default value is 0
%                 (i.e. recomputation disabled).
%                 
% OPTIONS.number_of_initial_iterations -
%                 Number of times the candidate regions are initially
%                 regrown. Default value is 2.
%                 
% OPTIONS.motion_models -
%                 Which motion model(s) to use. These are encoded as
%                 integers with:
%                 1 - constant motion
%                 2 - affine motion
%                 4 - eight parameter motion
%                 It is also possible to add the values for two or all
%                 three models. Then multiple candidate regions are
%                 created at each starting point, one for each motion
%                 model. This extension is not discussed in the thesis.
%                 Default value is 2, i.e. only affine motion.
%
% OPTIONS.verbose -
%                 If non-zero, trace messages are written while the
%                 algorithm is running. Probably only useful for
%                 debugging. Default is 0.
%                 
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

error('VELOCITY_SEGMENT is implemented as a mex-file. It has not been compiled on this platform.')
