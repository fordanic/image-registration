function result = conv3(signal, kernel, region_of_interest)
% CONV3
% 
% Perform convolution (actually correlation, see note 3 below) in up to
% 4 dimensions. The signal may be complex but the kernel must be real.
% 
% Call:
% result = conv3(signal, kernel)
% or
% result = conv3(signal, kernel, region_of_interest)
% 
% signal             - signal values
% kernel             - filter kernel
% region_of_interest - where to compute the convolution
% result             - filter output
% 
% Formats:
% signal is an up to 4-dimensional array.
% kernel is an up to 4-dimensional array.
% region of interest is an Nx2 matrix, where each row gives start
%	  and end indices for each signal dimension. N must be the
%	  same as the number of non-trailing singleton dimensions
%	  for signal.
% result is an up to 4-dimensional array. The size is the same as
%	  for the signal, unless region_of_interest is specified.
%	  The signal is assumed to be surrounded by zeros.
% 
% Note 1: Only double, nonsparse arrays are currently supported.
% Note 2: If signal or kernel has fewer than 3 dimensions,
%         trailing singleton dimensions are added.
% Note 3: Actually, the name of this function is misleading since it
%         doesn't do a proper convolution. The kernel is not mirrored.
% Note 4: The name is also misleading in that it can handle up to 4
%         dimensions. This is for historical reasons.
% 
% 
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

error('CONV3 is implemented as a mex-file. It has not been compiled on this platform.')
