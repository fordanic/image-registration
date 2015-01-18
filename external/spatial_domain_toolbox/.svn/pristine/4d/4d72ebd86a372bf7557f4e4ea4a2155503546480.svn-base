function [A1,b1,c1,A2,b2,c2,params] = make_Abc_interlace(signal,spatial_size,region_of_interest,options)
% [A1, b1, c1, A2, b2, c2] = make_Abc_interlace(signal, spatial_size, ...
%                                               region_of_interest, options)
% 
% Compute A, b, and c parameters for an interlaced 2D image. The
% parameters relate to the local signal model
% f(x) = x^T A x + b^T x + c
% and are determined by a Gaussian weighted least squares fit. This
% implementation uses a fast hierarchical scheme of separable filters,
% described in chapter 8 of "Polynomial Expansion for Orientation and
% Motion Estimation" by Gunnar Farnebäck.
% 
% signal                        - Signal values. Must be real and nonsparse 
%                                 and the number of dimensions must be two.
%
% spatial_size [optional]       - Size of the spatial support of the filters
%                                 along each dimension. Default value is 9.
%
% region_of_interest [optional] - An Nx2 matrix where each row contains
%                                 start and stop indices along the
%                                 corresponding dimensions. Default value
%                                 is all of the signal.
%
% options [optional]            - Struct array that may contain various
%                                 parameters that affect the algorithm.
%                                 These are explained below.
%
% A1                            - Computed A matrices for the first field.
%                                 A1 has 4 dimensions, where the first 2
%                                 indices indicate the position in the
%                                 signal and the last two contains the
%                                 matrix for each point.
%
% b1                            - Computed b vectors for the first field.
%                                 b1 has 3 dimensions, where the first 2
%                                 indices indicate the position in the
%                                 signal and the last one contains the
%                                 vector for each point.
%
% c1                            - Computed c scalars for the first field.
%                                 c1 has 2 dimensions.
%
% A2, b2, c2                    - Corresponding for the second field.
%
% params                        - Struct array containing the parameters
%                                 that has been used by the algorithm.
%
%
% The following fields may be specified in the options parameter:
%   
% options.sigma - Standard deviation of a Gaussian applicability. The
%                 default value is 0.15(K-1), where K is the spatial_size.
%                 However, if options.delta is set, that value is used
%                 instead.
%
% options.delta - The value of the gaussian applicability when it reaches
%                 the end of the supporting region along an axis. If both
%                 options.sigma and options.delta are set, the former is
%                 used.
%
% options.c     - Certainty mask. Must be spatially invariant and symmetric
%                 with respect to all axes and have a size compatible with
%                 the signal dimension and the spatial_size parameter.
%                 Default value is all ones.
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

N = ndims(signal);

if nargin < 2
    spatial_size = 9;
end

if spatial_size < 1
    error('What use would such a small kernel be?')
elseif mod(spatial_size, 2) ~= 1
    spatial_size = 2*floor((spatial_size-1)/2)+1;
    warning(sprintf('Only kernels of odd size are allowed. Changed the size to %d.', spatial_size))
end

if nargin < 3
    if N == 1
	region_of_interest = [1 size(signal, 1)];
    else
	region_of_interest = [ones(N, 1), size(signal)'];
    end
end

sigma = 0.15 * (spatial_size - 1);
certainty = ones([repmat(spatial_size, [1 N]) 1]);

if nargin == 4
    if isfield(options, 'sigma')
	sigma = options.sigma;
    elseif isfield(options, 'delta')
	sigma = n/sqrt(-2*log(delta));
    end
    if isfield(options, 'c')
	certainty = options.c;
    end
end

n = (spatial_size-1)/2;
a = exp(-(-n:n).^2/(2*sigma^2))';

switch N
    case 2
	% Set up applicability and basis functions.
	applicability = a*a';
	[x, y] = ndgrid(-n:n);
	b = cat(3, ones(size(x)), x, y, x.*x, y.*y, x.*y);
	nb = size(b, 3);

	% Compute the inverse metrics.
	c_center = zeros(spatial_size, spatial_size);
	c_center(1+mod(n,2):2:end, :) = 1;
	c_between = 1 - c_center;
	
	Q_center = zeros(nb,nb);
	for i = 1:nb
	    for j = i:nb
		Q_center(i,j) = sum(sum(b(:,:,i).*applicability.*c_center.*b(:,:,j)));
		Q_center(j,i) = Q_center(i,j);
		Q_between(i,j) = sum(sum(b(:,:,i).*applicability.*c_between.*b(:,:,j)));
		Q_between(j,i) = Q_between(i,j);
	    end
	end
	clear b applicability x y
	Qinv_center = inv(Q_center);
	Qinv_between = inv(Q_between);
	
	% Convolutions in the y-direction.
	kernely0 = a';
	kernely1 = (-n:n).*a';
	kernely2 = (-n:n).^2.*a';
	roiy = region_of_interest+[-n n;0 0];
	roiy(:,1) = max(roiy(:,1),ones(2,1));
	roiy(:,2) = min(roiy(:,2),size(signal)');
	conv_y0 = conv3(signal,kernely0,roiy);
	conv_y1 = conv3(signal,kernely1,roiy);
	conv_y2 = conv3(signal,kernely2,roiy);
	
	% Convolutions in the x-direction.
	kernelx0_center  = kernely0(:) .* c_center(:,1);
	kernelx1_center  = kernely1(:) .* c_center(:,1);
	kernelx2_center  = kernely2(:) .* c_center(:,1);
	kernelx0_between = kernely0(:) .* c_between(:,1);
	kernelx1_between = kernely1(:) .* c_between(:,1);
	kernelx2_between = kernely2(:) .* c_between(:,1);
	roix = region_of_interest;
	roix = roix(1:ndims(conv_y0),:);
	roix(2:end,:) = roix(2:end,:)+1-repmat(roix(2:end,1),[1 2]);
	conv_results_center = zeros([diff(region_of_interest')+1 6]);
	conv_results_between = zeros([diff(region_of_interest')+1 6]);
	conv_results_center(:,:,1) = conv3(conv_y0,kernelx0_center,roix); % y0x0
	conv_results_center(:,:,2) = conv3(conv_y0,kernelx1_center,roix); % y0x1
	conv_results_center(:,:,4) = conv3(conv_y0,kernelx2_center,roix); % y0x2
	conv_results_between(:,:,1) = conv3(conv_y0,kernelx0_between,roix); % y0x0
	conv_results_between(:,:,2) = conv3(conv_y0,kernelx1_between,roix); % y0x1
	conv_results_between(:,:,4) = conv3(conv_y0,kernelx2_between,roix); % y0x2
	clear conv_y0
	conv_results_center(:,:,3) = conv3(conv_y1,kernelx0_center,roix); % y1x0
	conv_results_center(:,:,6) = conv3(conv_y1,kernelx1_center,roix)/2; % y1x1
	conv_results_between(:,:,3) = conv3(conv_y1,kernelx0_between,roix); % y1x0
	conv_results_between(:,:,6) = conv3(conv_y1,kernelx1_between,roix)/2; % y1x1
	clear conv_y1
	conv_results_center(:,:,5) = conv3(conv_y2,kernelx0_center,roix); % y2x0
	conv_results_between(:,:,5) = conv3(conv_y2,kernelx0_between,roix); % y2x0
	clear conv_y2
	
	% Apply the inverse metrics.
	tmp = Qinv_center(1,1)*conv_results_center(:,:,1) + ...
	      Qinv_center(1,4)*conv_results_center(:,:,4) + ...
	      Qinv_center(1,5)*conv_results_center(:,:,5);
	conv_results_center(:,:,2) = Qinv_center(2,2)*conv_results_center(:,:,2);
	conv_results_center(:,:,3) = Qinv_center(3,3)*conv_results_center(:,:,3);
	conv_results_center(:,:,4) = Qinv_center(4,4)*conv_results_center(:,:,4) + ...
			      Qinv_center(4,1)*conv_results_center(:,:,1);
	conv_results_center(:,:,5) = Qinv_center(5,5)*conv_results_center(:,:,5) + ...
			      Qinv_center(5,1)*conv_results_center(:,:,1);
	conv_results_center(:,:,6) = Qinv_center(6,6)*conv_results_center(:,:,6);
	conv_results_center(:,:,1) = tmp;
	
	tmp = Qinv_between(1,1)*conv_results_between(:,:,1) + ...
	      Qinv_between(1,4)*conv_results_between(:,:,4) + ...
	      Qinv_between(1,5)*conv_results_between(:,:,5);
	conv_results_between(:,:,2) = Qinv_between(2,2)*conv_results_between(:,:,2);
	conv_results_between(:,:,3) = Qinv_between(3,3)*conv_results_between(:,:,3);
	conv_results_between(:,:,4) = Qinv_between(4,4)*conv_results_between(:,:,4) + ...
			      Qinv_between(4,1)*conv_results_between(:,:,1);
	conv_results_between(:,:,5) = Qinv_between(5,5)*conv_results_between(:,:,5) + ...
			      Qinv_between(5,1)*conv_results_between(:,:,1);
	conv_results_between(:,:,6) = Qinv_between(6,6)*conv_results_between(:,:,6);
	conv_results_between(:,:,1) = tmp;
	
	% Build A, b, and c.
	%
	% A, b, and c are obtained from the convolution results according to
	% 
	%   [4  6]    [2]
	% A=[6  5], b=[3], c=[1].
	%

	A1 = zeros([diff(region_of_interest')+1 2 2]);
	b1 = zeros([diff(region_of_interest')+1 2]);
	c1 = zeros(diff(region_of_interest')+1);
	A2 = zeros([diff(region_of_interest')+1 2 2]);
	b2 = zeros([diff(region_of_interest')+1 2]);
	c2 = zeros(diff(region_of_interest')+1);
	
	A1(1:2:end,:,1,1) = conv_results_center(1:2:end,:,4);
	A1(2:2:end,:,1,1) = conv_results_between(2:2:end,:,4);
	A2(1:2:end,:,1,1) = conv_results_between(1:2:end,:,4);
	A2(2:2:end,:,1,1) = conv_results_center(2:2:end,:,4);

	A1(1:2:end,:,2,2) = conv_results_center(1:2:end,:,5);
	A1(2:2:end,:,2,2) = conv_results_between(2:2:end,:,5);
	A2(1:2:end,:,2,2) = conv_results_between(1:2:end,:,5);
	A2(2:2:end,:,2,2) = conv_results_center(2:2:end,:,5);

	A1(1:2:end,:,1,2) = conv_results_center(1:2:end,:,6);
	A1(2:2:end,:,1,2) = conv_results_between(2:2:end,:,6);
	A2(1:2:end,:,1,2) = conv_results_between(1:2:end,:,6);
	A2(2:2:end,:,1,2) = conv_results_center(2:2:end,:,6);

	A1(:,:,2,1) = A1(:,:,1,2);
	A2(:,:,2,1) = A2(:,:,1,2);
	
	b1(1:2:end,:,1) = conv_results_center(1:2:end,:,2);
	b1(2:2:end,:,1) = conv_results_between(2:2:end,:,2);
	b2(1:2:end,:,1) = conv_results_between(1:2:end,:,2);
	b2(2:2:end,:,1) = conv_results_center(2:2:end,:,2);

	b1(1:2:end,:,2) = conv_results_center(1:2:end,:,3);
	b1(2:2:end,:,2) = conv_results_between(2:2:end,:,3);
	b2(1:2:end,:,2) = conv_results_between(1:2:end,:,3);
	b2(2:2:end,:,2) = conv_results_center(2:2:end,:,3);

	c1(1:2:end,:) = conv_results_center(1:2:end,:,1);
	c1(2:2:end,:) = conv_results_between(2:2:end,:,1);
	c2(1:2:end,:) = conv_results_between(1:2:end,:,1);
	c2(2:2:end,:) = conv_results_center(2:2:end,:,1);

	A1 = squeeze(A1);
	A2 = squeeze(A2);
	b1 = squeeze(b1);
	b2 = squeeze(b2);
	c1 = squeeze(c1);
	c2 = squeeze(c2);
	
	if mod(region_of_interest(1,1), 2) == 0
	   tmp = A1;
	   A1 = A2;
	   A2 = tmp;
	   tmp = b1;
	   b1 = b2;
	   b2 = tmp;
	   tmp = c1;
	   c1 = c2;
	   c2 = tmp;
	end
	
    otherwise
	error('Only 2D signals supported by this function.')
end

if nargout>6
    params.spatial_size = spatial_size;
    params.region_of_interest = region_of_interest;
    params.sigma = sigma;
    params.delta = delta;
end
