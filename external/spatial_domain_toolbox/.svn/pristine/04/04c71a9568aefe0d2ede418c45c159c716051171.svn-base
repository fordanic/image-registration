function [r, cout] = normconv(signal, certainty, basis, applicability, ...
			      region_of_interest, options)
% R = NORMCONV(SIGNAL, CERTAINTY, BASIS, APPLICABILITY, ...
%              REGION_OF_INTEREST, OPTIONS)
% or
% [R COUT] = NORMCONV(SIGNAL, CERTAINTY, BASIS, APPLICABILITY, ...
%                     REGION_OF_INTEREST, OPTIONS)
% 
% Perform normalized convolution. Complex signals and basis functions are
% supported, at arbitrary dimensionality. The method is described in chapter
% 3 of Gunnar Farnebäck's thesis "Polynomial Expansion for Orientation and
% Motion Estimation".
% 
% SIGNAL
%   Signal values. The number of dimensions, N, is arbitrary.
%
% CERTAINTY
%   Certainty values. Must be real with size identical to that of the signal.
%
% BASIS
%   Subspace basis functions. Each basis function must have the same
%   dimensionality as the signal. Then they are stacked along dimension N+1.
%
% APPLICABILITY
%   Applicability function for the basis. Must be real with size
%   identical to that of each separate basis function.
%
% REGION_OF_INTEREST [optional]
%   Specification of where to compute results. This can have two formats.
%   1. An Nx2 matrix where each row contains start and stop
%      indices along the corresponding dimensions (box region of interest).
%   2. An array of the same size as the signal, where non-zero values
%      indicate that the results should be computed (mask region of
%      interest).
%   Default value is all of the signal.
%
% OPTIONS [optional]
%   Struct array that may contain various parameters that affect the
%   algorithm. These are explained below. The default value is no options.
%
% R
%   Computed coordinates for the basis functions. R has N+1 dimensions,
%   where the first N indices indicate the position in the signal and the
%   last dimension holds the coordinates. In the case that
%   REGION_OF_INTEREST is of the box type and less than N-dimensional,
%   the singleton dimensions are removed.
% 
% COUT
%   Output certainty. This is only available if OPTIONS includes a
%   field cout_func. COUT has N+K dimensions, where the first N
%   indices indicate the position in the signal and the last K dimensions
%   hold the output certainty. In the case that REGION_OF_INTEREST
%   is of the box type and less than N-dimensional, the singleton
%   dimensions are removed.
%
%
% The optional parameters can be omitted freely since it is always possible
% to distinguish them.
% 
% The following fields may be specified in the OPTIONS parameter:
% OPTIONS.region_of_interest -
%                 Same as the REGION_OF_INTEREST parameter. Both can not
%                 be used in the same call.
%
% OPTIONS.non_optimized -
%                 Do not use any optimizations if this field exists
%                 and is non-zero.
%
% OPTIONS.cout_func -
%                 Name of a function to compute output certainty.
%                 The function will be called for each point where
%                 coordinates are computed. The call is
%                 of the form cout_func(G, G0, h, r, cout_data)
%                 where G and G0 are as in equation (3.18),
%                 r is a vector with the coordinates,
%                 h equals G*r, and cout_data is specified below.
%                 The output from cout_func must be a numeric array
%                 of the same size at all points.
%                 
% OPTIONS.cout_data -
%                 Arbitrary data passed on to cout_func.
%
%
% Note 1: Only double, nonsparse arrays are currently supported.
%
% Note 2: Trailing singleton dimensions cannot exist in Matlab.
%         If there is a mismatch in the number of dimensions of the
%         parameters it will be assumed that there are additional
%         singleton dimensions, if that turns out to make sense.
%         Particularly, in the case of a single basis function, basis
%         will normally be N-dimensional instead of (N+1)-dimensional.
%         The same goes for result.
% 
%         In the case of column vectors, the second singleton dimension
%         will be ignored, although it is included by Matlab.
%
% Note 3: The special cases of 1-3 dimensions have been optimized for
%         speed. To use the general algorithm also in these cases,
%         set the OPTIONS.non_optimized field. This will certainly be slower
%         but can be used to verify that the optimized code works as
%         intended. The results should not differ by more than
%         numerical deviations.
%
% Note 4: For the box region of interest the spatial dimensions of the
%         output have the size given by the region of interest. For the
%         mask region of interest the spatial dimensions of the output
%         are the same as for the signal. Non-computed output values
%         are set to zero.
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden

error('NORMCONV is implemented as a mex-file. It has not been compiled on this platform.')
