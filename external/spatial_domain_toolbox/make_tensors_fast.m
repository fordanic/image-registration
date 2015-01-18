function [T, params] = make_tensors_fast(signal, spatial_size, ...
					 region_of_interest, options)
% T = MAKE_TENSORS_FAST(SIGNAL, SPATIAL_SIZE, REGION_OF_INTEREST, OPTIONS)
% or
% [T,PARAMS]=MAKE_TENSORS_FAST(SIGNAL,SPATIAL_SIZE,REGION_OF_INTEREST,OPTIONS)
% 
% Compute orientation tensors in up to four dimensions. The tensors are
% computed according to an algorithm described in chapter 5 of Gunnar
% Farnebäck's thesis, "Polynomial Expansion for Orientation and Motion
% Estimation". This implementation uses the "Separable Convolution" method
% with completely separable filters in a hierarchical scheme. The
% applicability (weighting function) is limited to be a Gaussian.
% 
% SIGNAL                        - Signal values. Must be real and nonsparse 
%                                 and the number of dimensions, N, must be
%                                 at most four.
%
% SPATIAL_SIZE [optional]       - Size of the spatial support of the filters
%                                 along each dimension. Default value is 9.
%
% REGION_OF_INTEREST [optional] - An Nx2 matrix where each row contains
%                                 start and stop indices along the
%                                 corresponding dimensions. Default value
%                                 is all of the signal.
%
% OPTIONS [optional]            - Struct array that may contain various
%                                 parameters that affect the algorithm.
%                                 These are explained below.
%
% T                             - Computed tensor field. T has N+2
%                                 dimensions, where the first N indices
%                                 indicate the position in the signal and
%                                 the last two contain the tensor for each
%                                 point. In the case that REGION_OF_INTEREST
%                                 is less than N-dimensional, the singleton
%                                 dimensions are removed.
%
% PARAMS                        - Struct array containing the parameters
%                                 that have been used by the algorithm.
%
%
% The following fields may be specified in the OPTIONS parameter:
% OPTIONS.gamma - Relation between the contribution to the tensor from the
%                 linear and the quadratic parts of the signal, as specified
%                 in equation (5.19). 0 means that only the quadratic part
%                 matters while a very large number means that only the
%                 linear part is used. Default value is 1/(8*sigma^2).
%
% OPTIONS.sigma - Standard deviation of a Gaussian applicability. The
%                 default value is 0.15(K-1), where K is the SPATIAL_SIZE.
%
% OPTIONS.delta - The value of the gaussian applicability when it reaches
%                 the end of the supporting region along an axis. If both
%                 OPTIONS.sigma and OPTIONS.delta are set, the former is
%                 used.
%
% OPTIONS.c     - Certainty mask. Must be spatially invariant and symmetric
%                 with respect to all axes and have a size compatible with
%                 the signal dimension and the SPATIAL_SIZE parameter.
%                 Default value is all ones. One application of this option
%                 is in conjunction with interlaced images.
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

% We are going to modify the value returned by nargin, so we copy it to a
% variable.
numin = nargin;

N = ndims(signal);
if N == 2 & size(signal, 2) == 1
    N = 1;
end

if numin < 2 | (numin == 2 & isstruct(spatial_size))
    if numin == 2
	numin = 3;
	region_of_interest = spatial_size;
    end
    spatial_size = 9;
end
    
if spatial_size < 1
    error('What use would such a small kernel be?')
elseif mod(spatial_size, 2)~=1
    spatial_size = 2*floor((spatial_size-1)/2) + 1;
    warning(sprintf('Only kernels of odd size are allowed. Changed the size to %d.', spatial_size))
end

if numin < 3 | (numin == 3 & isstruct(region_of_interest))
    if numin == 3
	numin = 4;
	options = region_of_interest;
    end
    if N == 1
	region_of_interest = [1 size(signal, 1)];
    else
	region_of_interest = [ones(N, 1), size(signal)'];
    end
end

sigma = 0.15 * (spatial_size - 1);
gamma = -1; % Mark gamma as uninitialized.
certainty = ones(repmat(spatial_size, [1 N]));

if numin == 4
    if isfield(options, 'sigma')
	sigma = options.sigma;
    elseif isfield(options, 'delta')
	sigma = n/sqrt(-2*log(delta));
    end
    if isfield(options, 'gamma')
	gamma = options.gamma;
    end
    if isfield(options, 'c')
	certainty = options.c;
    end
end

n = (spatial_size - 1) / 2;
a = exp(-(-n:n).^2/(2*sigma^2))';

% If gamma has not been set explicitly, use default value.
if gamma < 0
    gamma = 1 / (8*sigma^2);
end

switch N
    case 1
      	% Comment: Orientation tensors in 1D are fairly pointless and only
      	%          included here for completeness.
      	%
	% Set up applicability and basis functions.
	applicability = a;
	x = (-n:n)';
	b = [ones(size(x)) x x.*x];
	nb = size(b, 2);

	% Compute the inverse metric.
	Q = zeros(nb, nb);
	for i = 1:nb
	    for j = i:nb
		Q(i,j) = sum(b(:,i).*applicability.*certainty.*b(:,j));
		Q(j,i) = Q(i,j);
	    end
	end
	clear b applicability x y
	Qinv = inv(Q);
	
	% Convolutions in the x-direction.
	kernelx0 = a;
	kernelx1 = (-n:n)'.*a;
	kernelx2 = (-n:n)'.^2.*a;
	roix = region_of_interest;
	roix(1) = max(roix(1), 1);
	roix(2) = min(roix(2), length(signal));
	conv_results = zeros([diff(region_of_interest')+1 3]);
	conv_results(:,1) = conv3(signal, kernelx0, roix);
	conv_results(:,2) = conv3(signal, kernelx1, roix);
	conv_results(:,3) = conv3(signal, kernelx2, roix);
	
	% Apply the inverse metric.
	conv_results(:,2) = Qinv(2,2) * conv_results(:,2);
	conv_results(:,3) = Qinv(3,3) * conv_results(:,3) + ...
			    Qinv(3,1) * conv_results(:,1);
	
	% Build tensor components.
	%
	% It's more efficient in matlab code to do a small matrix
	% multiplication "manually" in parallell over all the points
	% than doing a multiple loop over the points and computing the
	% matrix products "automatically".
	% The tensor is of the form A*A'+gamma*b*b', where A and b are
	% composed from the convolution results according to
	% 
	% A=[3], b=[2].
	%
	% Thus (excluding gamma)
	%
	% T=[3*3+2*2].

	T = zeros([diff(region_of_interest')+1 1 1]);
	T(:,1,1) = conv_results(:,3).^2 + gamma * conv_results(:,2).^2;
	T = squeeze(T);
	
    
    case 2
	% Set up applicability and basis functions.
	applicability = a*a';
	[x,y] = ndgrid(-n:n);
	b = cat(3, ones(size(x)), x, y, x.*x, y.*y, x.*y);
	nb = size(b, 3);

	% Compute the inverse metric.
	Q = zeros(nb, nb);
	for i = 1:nb
	    for j = i:nb
		Q(i,j) = sum(sum(b(:,:,i).*applicability.*certainty.*b(:,:,j)));
		Q(j,i) = Q(i,j);
	    end
	end
	clear b applicability x y
	Qinv = inv(Q);
	
	% Convolutions in the y-direction.
	kernely0 = a';
	kernely1 = (-n:n).*a';
	kernely2 = (-n:n).^2.*a';
	roiy = region_of_interest + [-n n;0 0];
	roiy(:,1) = max(roiy(:,1), ones(2,1));
	roiy(:,2) = min(roiy(:,2), size(signal)');
	conv_y0 = conv3(signal, kernely0, roiy);
	conv_y1 = conv3(signal, kernely1, roiy);
	conv_y2 = conv3(signal, kernely2, roiy);
	
	% Convolutions in the x-direction.
	kernelx0 = kernely0(:);
	kernelx1 = kernely1(:);
	kernelx2 = kernely2(:);
	roix = region_of_interest;
	roix = roix(1:ndims(conv_y0),:);
	roix(2:end,:) = roix(2:end,:) + 1 - repmat(roix(2:end,1), [1 2]);
	conv_results = zeros([diff(region_of_interest')+1 6]);
	conv_results(:,:,1) = conv3(conv_y0, kernelx0, roix); % y0x0
	conv_results(:,:,2) = conv3(conv_y0, kernelx1, roix); % y0x1
	conv_results(:,:,4) = conv3(conv_y0, kernelx2, roix); % y0x2
	clear conv_y0
	conv_results(:,:,3) = conv3(conv_y1, kernelx0, roix); % y1x0
	conv_results(:,:,6) = conv3(conv_y1, kernelx1, roix)/2; % y1x1
	clear conv_y1
	conv_results(:,:,5) = conv3(conv_y2, kernelx0, roix); % y2x0
	clear conv_y2
	
	% Apply the inverse metric.
	conv_results(:,:,2) = Qinv(2,2) * conv_results(:,:,2);
	conv_results(:,:,3) = Qinv(3,3) * conv_results(:,:,3);
	conv_results(:,:,4) = Qinv(4,4) * conv_results(:,:,4) + ...
			      Qinv(4,1) * conv_results(:,:,1);
	conv_results(:,:,5) = Qinv(5,5) * conv_results(:,:,5) + ...
			      Qinv(5,1) * conv_results(:,:,1);
	conv_results(:,:,6) = Qinv(6,6) * conv_results(:,:,6);
	
	% Build tensor components.
	%
	% It's more efficient in matlab code to do a small matrix
	% multiplication "manually" in parallell over all the points
	% than doing a multiple loop over the points and computing the
	% matrix products "automatically".
	% The tensor is of the form A*A'+gamma*b*b', where A and b are
	% composed from the convolution results according to
	% 
	%   [4  6]    [2]
	% A=[6  5], b=[3].
	%
	% Thus (excluding gamma)
	%
	%   [4*4+6*6+2*2 4*6+5*6+2*3]
	% T=[4*6+5*6+2*3 6*6+5*5+3*3].

	T = zeros([diff(region_of_interest')+1 2 2]);
	T(:,:,1,1) = conv_results(:,:,4).^2+...
		     conv_results(:,:,6).^2+...
		     gamma * conv_results(:,:,2).^2;
	T(:,:,2,2) = conv_results(:,:,6).^2+...
		     conv_results(:,:,5).^2+...
		     gamma * conv_results(:,:,3).^2;
	T(:,:,1,2) = (conv_results(:,:,4)+...
		      conv_results(:,:,5)).*...
		     conv_results(:,:,6)+...
		     gamma * conv_results(:,:,2).*...
		     conv_results(:,:,3);
	T(:,:,2,1) = T(:,:,1,2);
	
	T = squeeze(T);
	
    case 3
	% Set up applicability and basis functions.
	applicability = outerprod(a, a, a);
	[x,y,t] = ndgrid(-n:n);
	b = cat(4, ones(size(x)), x, y, t, x.*x, y.*y, t.*t, x.*y, x.*t, y.*t);
	nb = size(b, 4);

	% Compute the inverse metric.
	Q = zeros(nb, nb);
	for i = 1:nb
	    for j = i:nb
		Q(i,j) = sum(sum(sum(b(:,:,:,i).*applicability.*certainty.*b(:,:,:,j))));
		Q(j,i) = Q(i,j);
	    end
	end
	clear b applicability x y t
	Qinv = inv(Q);

	% Convolutions in the t-direction
	kernelt0 = reshape(a,[1 1 spatial_size]);
	kernelt1 = reshape((-n:n)'.*a,[1 1 spatial_size]);
	kernelt2 = reshape(((-n:n).^2)'.*a,[1 1 spatial_size]);
	roit = region_of_interest+[-n n;-n n;0 0];
	roit(:,1) = max(roit(:,1), ones(3,1));
	roit(:,2) = min(roit(:,2), size(signal)');
	conv_t0 = conv3(signal, kernelt0, roit);
	conv_t1 = conv3(signal, kernelt1, roit);
	conv_t2 = conv3(signal, kernelt2, roit);
	
	% Convolutions in the y-direction
	kernely0 = reshape(kernelt0, [1 spatial_size]);
	kernely1 = reshape(kernelt1, [1 spatial_size]);
	kernely2 = reshape(kernelt2, [1 spatial_size]);
	roiy = region_of_interest + [-n n;0 0;0 0];
	roiy(:,1) = max(roiy(:,1), ones(3,1));
	roiy(:,2) = min(roiy(:,2), size(signal)');
	if diff(roiy(3,:)) == 0
	    roiy = roiy(1:2,:);
	else
	    roiy(3,:) = roiy(3,:) + 1 - roiy(3,1);
	end
	conv_t0y0 = conv3(conv_t0, kernely0, roiy);
	conv_t0y1 = conv3(conv_t0, kernely1, roiy);
	conv_t0y2 = conv3(conv_t0, kernely2, roiy);
	clear conv_t0
	conv_t1y0 = conv3(conv_t1, kernely0, roiy);
	conv_t1y1 = conv3(conv_t1, kernely1, roiy);
	clear conv_t1
	conv_t2y0 = conv3(conv_t2, kernely0, roiy);
	clear conv_t2
	
	% Convolutions in the x-direction
	kernelx0 = reshape(kernelt0, [spatial_size 1]);
	kernelx1 = reshape(kernelt1, [spatial_size 1]);
	kernelx2 = reshape(kernelt2, [spatial_size 1]);
	roix = region_of_interest;
	roix = roix(1:ndims(conv_t0y0),:);
	roix(2:end,:) = roix(2:end,:) + 1 - repmat(roix(2:end,1),[1 2]);
	conv_results = zeros([diff(region_of_interest')+1 10]);
	conv_results(:,:,:,1) = conv3(conv_t0y0, kernelx0, roix); % t0y0x0
	conv_results(:,:,:,2) = conv3(conv_t0y0, kernelx1, roix); % t0y0x1
	conv_results(:,:,:,5) = conv3(conv_t0y0, kernelx2, roix); % t0y0x2
	clear conv_t0y0
	conv_results(:,:,:,3) = conv3(conv_t0y1, kernelx0, roix); % t0y1x0
	conv_results(:,:,:,8) = conv3(conv_t0y1, kernelx1, roix)/2; % t0y1x1
	clear conv_t0y1
	conv_results(:,:,:,6) = conv3(conv_t0y2, kernelx0, roix); % t0y2x0
	clear conv_t0y2
	conv_results(:,:,:,4) = conv3(conv_t1y0, kernelx0, roix); % t1y0x0
	conv_results(:,:,:,9) = conv3(conv_t1y0, kernelx1, roix)/2; % t1y0x1
	clear conv_t1y0
	conv_results(:,:,:,10) = conv3(conv_t1y1, kernelx0, roix)/2; % t1y1x0
	clear conv_t1y1
	conv_results(:,:,:,7) = conv3(conv_t2y0, kernelx0, roix); % t2y0x0
	clear conv_t2y0
	
	% Apply the inverse metric.
	conv_results(:,:,:,2) = Qinv(2,2) * conv_results(:,:,:,2);
	conv_results(:,:,:,3) = Qinv(3,3) * conv_results(:,:,:,3);
	conv_results(:,:,:,4) = Qinv(4,4) * conv_results(:,:,:,4);
	conv_results(:,:,:,5) = Qinv(5,5) * conv_results(:,:,:,5) + ...
				Qinv(5,1) * conv_results(:,:,:,1);
	conv_results(:,:,:,6) = Qinv(6,6) * conv_results(:,:,:,6) + ...
				Qinv(6,1) * conv_results(:,:,:,1);
	conv_results(:,:,:,7) = Qinv(7,7) * conv_results(:,:,:,7) + ...
				Qinv(7,1) * conv_results(:,:,:,1);
	conv_results(:,:,:,8) = Qinv(8,8) * conv_results(:,:,:,8);
	conv_results(:,:,:,9) = Qinv(9,9) * conv_results(:,:,:,9);
	conv_results(:,:,:,10) = Qinv(10,10) * conv_results(:,:,:,10);

	% Build tensor components.
	%
	% It's more efficient in matlab code to do a small matrix
	% multiplication "manually" in parallell over all the points
	% than doing a multiple loop over the points and computing the
	% matrix products "automatically".
	% The tensor is of the form T=A*A'+gamma*b*b', where A and b are
	% composed from the convolution results according to
	% 
	%   [5  8  9]    [2]
	% A=[8  6 10], b=[3].
	%   [9 10  7]    [4]
	%
	% Thus (excluding gamma)
	%
	%   [5*5+8*8+9*9+2*2  5*8+6*8+9*10+2*3  5*9+8*10+7*9+2*4 ]
	% T=[5*8+6*8+9*10+2*3 8*8+6*6+10*10+3*3 8*9+6*10+7*10+3*4].
	%   [5*9+8*10+7*9+2*4 8*9+6*10+7*10+3*4 9*9+10*10+7*7+4*4]
	
	T = zeros([diff(region_of_interest')+1 3 3]);
	T(:,:,:,1,1) = conv_results(:,:,:,5).^2+...
		       conv_results(:,:,:,8).^2+...
		       conv_results(:,:,:,9).^2+...
		       gamma * conv_results(:,:,:,2).^2;
	T(:,:,:,2,2) = conv_results(:,:,:,8).^2+...
		       conv_results(:,:,:,6).^2+...
		       conv_results(:,:,:,10).^2+...
		       gamma*conv_results(:,:,:,3).^2;
	T(:,:,:,3,3) = conv_results(:,:,:,9).^2+...
		       conv_results(:,:,:,10).^2+...
		       conv_results(:,:,:,7).^2+...
		       gamma*conv_results(:,:,:,4).^2;
	T(:,:,:,1,2) = (conv_results(:,:,:,5)+...
			conv_results(:,:,:,6)).*...
		       conv_results(:,:,:,8)+...
		       conv_results(:,:,:,9).*...
		       conv_results(:,:,:,10)+...
		       gamma*conv_results(:,:,:,2).*...
		       conv_results(:,:,:,3);
	T(:,:,:,2,1) = T(:,:,:,1,2);
	T(:,:,:,1,3) = (conv_results(:,:,:,5)+...
			conv_results(:,:,:,7)).*...
		       conv_results(:,:,:,9)+...
		       conv_results(:,:,:,8).*...
		       conv_results(:,:,:,10)+...
		       gamma*conv_results(:,:,:,2).*...
		       conv_results(:,:,:,4);
	T(:,:,:,3,1) = T(:,:,:,1,3);
	T(:,:,:,2,3) = (conv_results(:,:,:,6)+...
			conv_results(:,:,:,7)).*...
		       conv_results(:,:,:,10)+...
		       conv_results(:,:,:,8).*...
		       conv_results(:,:,:,9)+...
		       gamma*conv_results(:,:,:,3).*...
		       conv_results(:,:,:,4);
	T(:,:,:,3,2) = T(:,:,:,2,3);
	
	T = squeeze(T);
	
    case 4
	% Set up applicability and basis functions.
	applicability = outerprod(a, a, a, a);
	[x,y,z,t] = ndgrid(-n:n);
	b = cat(5, ones(size(x)), x, y, z, t, x.*x, y.*y, z.*z, t.*t,...
		x.*y, x.*z, x.*t, y.*z, y.*t, z.*t);
	nb = size(b, 5);

	% Compute the inverse metric.
	Q = zeros(nb, nb);
	for i = 1:nb
	    for j = i:nb
		Q(i,j) = sum(sum(sum(sum(b(:,:,:,:,i).*applicability.*...
					 certainty.*b(:,:,:,:,j)))));
		Q(j,i) = Q(i,j);
	    end
	end
	clear b applicability x y z t
	Qinv = inv(Q);

	% Convolutions in the t-direction
	kernelt0 = reshape(a,[1 1 1  spatial_size]);
	kernelt1 = reshape((-n:n)'.*a, [1 1 1 spatial_size]);
	kernelt2 = reshape(((-n:n).^2)'.*a, [1 1 1 spatial_size]);
	roit = region_of_interest+[-n n;-n n;-n n;0 0];
	roit(:,1) = max(roit(:,1), ones(4,1));
	roit(:,2) = min(roit(:,2), size(signal)');
	conv_t0 = conv3(signal, kernelt0, roit);
	conv_t1 = conv3(signal, kernelt1, roit);
	conv_t2 = conv3(signal, kernelt2, roit);
	
	% Convolutions in the z-direction
	kernelz0 = reshape(kernelt0, [1 1 spatial_size]);
	kernelz1 = reshape(kernelt1, [1 1 spatial_size]);
	kernelz2 = reshape(kernelt2, [1 1 spatial_size]);
	roiz = region_of_interest+[-n n;-n n;0 0;0 0];
	roiz(:,1) = max(roiz(:,1), ones(4,1));
	roiz(:,2) = min(roiz(:,2), size(signal)');
	if diff(roiz(4,:)) == 0
	    roiz = roiz(1:2,:);
	else
	    roiz(4,:) = roiz(4,:) + 1 - roiz(4,1);
	end
	conv_t0z0 = conv3(conv_t0, kernelz0, roiz);
	conv_t0z1 = conv3(conv_t0, kernelz1, roiz);
	conv_t0z2 = conv3(conv_t0, kernelz2, roiz);
	clear conv_t0
	conv_t1z0 = conv3(conv_t1, kernelz0, roiz);
	conv_t1z1 = conv3(conv_t1, kernelz1, roiz);
	clear conv_t1
	conv_t2z0 = conv3(conv_t2, kernelz0, roiz);
	clear conv_t2
	
	% Convolutions in the y-direction
	kernely0 = reshape(kernelt0, [1 spatial_size]);
	kernely1 = reshape(kernelt1, [1 spatial_size]);
	kernely2 = reshape(kernelt2, [1 spatial_size]);
	roiy = region_of_interest + [-n n;0 0;0 0;0 0];
	roiy(:,1) = max(roiy(:,1), ones(4,1));
	roiy(:,2) = min(roiy(:,2), size(signal)');
	roiy = roiy(1:ndims(conv_t0z0),:);
	roiy(3:end,:) = roiy(3:end,:) + 1 - repmat(roiy(3:end,1),[1 2]);
	conv_t0z0y0 = conv3(conv_t0z0, kernely0, roiy);
	conv_t0z0y1 = conv3(conv_t0z0, kernely1, roiy);
	conv_t0z0y2 = conv3(conv_t0z0, kernely2, roiy);
	clear conv_t0z0
	conv_t0z1y0 = conv3(conv_t0z1, kernely0, roiy);
	conv_t0z1y1 = conv3(conv_t0z1, kernely1, roiy);
	clear conv_t0z1
	conv_t0z2y0 = conv3(conv_t0z2, kernely0, roiy);
	clear conv_t0z2
	conv_t1z0y0 = conv3(conv_t1z0, kernely0, roiy);
	conv_t1z0y1 = conv3(conv_t1z0, kernely1, roiy);
	clear conv_t1z0
	conv_t1z1y0 = conv3(conv_t1z1, kernely0, roiy);
	clear conv_t1z1
	conv_t2z0y0 = conv3(conv_t2z0, kernely0, roiy);
	clear conv_t2z0
	
	% Convolutions in the x-direction
	kernelx0 = reshape(kernelt0, [spatial_size 1]);
	kernelx1 = reshape(kernelt1, [spatial_size 1]);
	kernelx2 = reshape(kernelt2, [spatial_size 1]);
	roix = region_of_interest;
	roix = roix(1:ndims(conv_t0z0y0),:);
	roix(2:end,:) = roix(2:end,:) + 1 - repmat(roix(2:end,1),[1 2]);
	conv_results = zeros([diff(region_of_interest')+1 15]);
	conv_results(:,:,:,:,1) = conv3(conv_t0z0y0, kernelx0, roix); % t0z0y0x0
	conv_results(:,:,:,:,2) = conv3(conv_t0z0y0, kernelx1, roix); % t0z0y0x1
	conv_results(:,:,:,:,6) = conv3(conv_t0z0y0, kernelx2, roix); % t0z0y0x2
	clear conv_t0z0y0
	conv_results(:,:,:,:,3) = conv3(conv_t0z0y1, kernelx0, roix); % t0z0y1x0
	conv_results(:,:,:,:,10) = conv3(conv_t0z0y1, kernelx1, roix)/2; % t0z0y1x1
	clear conv_t0z0y1
	conv_results(:,:,:,:,7) = conv3(conv_t0z0y2, kernelx0, roix); % t0z0y2x0
	clear conv_t0z0y2
	conv_results(:,:,:,:,4) = conv3(conv_t0z1y0, kernelx0, roix); % t0z1y0x0
	conv_results(:,:,:,:,11) = conv3(conv_t0z1y0, kernelx1, roix)/2; % t0z1y0x1
	clear conv_t0z1y0
	conv_results(:,:,:,:,13) = conv3(conv_t0z1y1, kernelx0, roix)/2; % t0z1y1x0
	clear conv_t0z1y1
	conv_results(:,:,:,:,8) = conv3(conv_t0z2y0, kernelx0, roix); % t0z2y0x0
	clear conv_t0z2y0
	conv_results(:,:,:,:,5) = conv3(conv_t1z0y0, kernelx0, roix); % t1z0y0x0
	conv_results(:,:,:,:,12) = conv3(conv_t1z0y0, kernelx1, roix)/2; % t1z0y0x1
	clear conv_t1z0y0
	conv_results(:,:,:,:,14) = conv3(conv_t1z0y1, kernelx0, roix)/2; % t1z0y1x0
	clear conv_t1z0y1
	conv_results(:,:,:,:,15) = conv3(conv_t1z1y0, kernelx0, roix)/2; % t1z1y0x0
	clear conv_t1z1y0
	conv_results(:,:,:,:,9) = conv3(conv_t2z0y0, kernelx0, roix); % t2z0y0x0
	clear conv_t2z0y0
	
	% Apply the inverse metric.
	conv_results(:,:,:,:,2) = Qinv(2,2) * conv_results(:,:,:,:,2);
	conv_results(:,:,:,:,3) = Qinv(3,3) * conv_results(:,:,:,:,3);
	conv_results(:,:,:,:,4) = Qinv(4,4) * conv_results(:,:,:,:,4);
	conv_results(:,:,:,:,5) = Qinv(5,5) * conv_results(:,:,:,:,5);
	conv_results(:,:,:,:,6) = Qinv(6,6) * conv_results(:,:,:,:,6) + ...
				  Qinv(6,1) * conv_results(:,:,:,:,1);
	conv_results(:,:,:,:,7) = Qinv(7,7) * conv_results(:,:,:,:,7) + ...
				  Qinv(7,1) * conv_results(:,:,:,:,1);
	conv_results(:,:,:,:,8) = Qinv(8,8) * conv_results(:,:,:,:,8) + ...
				  Qinv(8,1) * conv_results(:,:,:,:,1);
	conv_results(:,:,:,:,9) = Qinv(9,9) * conv_results(:,:,:,:,9) + ...
				  Qinv(9,1) * conv_results(:,:,:,:,1);
	conv_results(:,:,:,:,10) = Qinv(10,10) * conv_results(:,:,:,:,10);
	conv_results(:,:,:,:,11) = Qinv(11,11) * conv_results(:,:,:,:,11);
	conv_results(:,:,:,:,12) = Qinv(12,12) * conv_results(:,:,:,:,12);
	conv_results(:,:,:,:,13) = Qinv(13,13) * conv_results(:,:,:,:,13);
	conv_results(:,:,:,:,14) = Qinv(14,14) * conv_results(:,:,:,:,14);
	conv_results(:,:,:,:,15) = Qinv(15,15) * conv_results(:,:,:,:,15);
	
	% Build tensor components.
	%
	% It's more efficient in matlab code to do a small matrix
	% multiplication "manually" in parallell over all the points
	% than doing a multiple loop over the points and computing the
	% matrix products "automatically".
	% The tensor is of the form T=A*A'+gamma*b*b', where A and b are
	% composed from the convolution results according to
	% 
	%   [6  10 11 12]    [2]
	%   [10  7 13 14]    [3] 
	% A=[11 13  8 15], b=[4].
	%   [12 14 15  9]    [5]
	%
	% Thus (excluding gamma)
	%
	%   [6*6+10*10+11*11+12*12+2*2 6*10+7*10+11*13+12*14+2*3
	%   [6*10+7*10+11*13+12*14+2*3 10*10+7*7+13*13+14*14+3*3
	% T=[6*11+10*13+8*11+12*15+2*4 10*11+7*13+8*13+14*15+3*4
	%   [6*12+10*14+11*15+9*12+2*5 10*12+7*14+13*15+9*14+3*5
	%
	%    6*11+10*13+8*11+12*15+2*4 6*12+10*14+11*15+9*12+2*5]
	%    10*11+7*13+8*13+14*15+3*4 10*12+7*14+13*15+9*14+3*5]
	%    11*11+13*13+8*8+15*15+4*4 11*12+13*14+8*15+9*15+4*5].
	%    11*12+13*14+8*15+9*15+4*5 12*12+14*14+15*15+9*9+5*5]
	
	T = zeros([diff(region_of_interest')+1 3 3]);
	T(:,:,:,:,1,1) = conv_results(:,:,:,:,6).^2+...
			 conv_results(:,:,:,:,10).^2+...
			 conv_results(:,:,:,:,11).^2+...
			 conv_results(:,:,:,:,12).^2+...
			 gamma * conv_results(:,:,:,:,2).^2;
	T(:,:,:,:,2,2) = conv_results(:,:,:,:,10).^2+...
			 conv_results(:,:,:,:,7).^2+...
			 conv_results(:,:,:,:,13).^2+...
			 conv_results(:,:,:,:,14).^2+...
			 gamma * conv_results(:,:,:,:,3).^2;
	T(:,:,:,:,3,3) = conv_results(:,:,:,:,11).^2+...
			 conv_results(:,:,:,:,13).^2+...
			 conv_results(:,:,:,:,8).^2+...
			 conv_results(:,:,:,:,15).^2+...
			 gamma * conv_results(:,:,:,:,4).^2;
	T(:,:,:,:,4,4) = conv_results(:,:,:,:,11).^2+...
			 conv_results(:,:,:,:,13).^2+...
			 conv_results(:,:,:,:,8).^2+...
			 conv_results(:,:,:,:,15).^2+...
			 gamma * conv_results(:,:,:,:,4).^2;
	T(:,:,:,:,1,2) = (conv_results(:,:,:,:,6)+...
			  conv_results(:,:,:,:,7)).*...
			 conv_results(:,:,:,:,10)+...
			 conv_results(:,:,:,:,11).*...
			 conv_results(:,:,:,:,13)+...
			 conv_results(:,:,:,:,12).*...
			 conv_results(:,:,:,:,14)+...
			 gamma * conv_results(:,:,:,:,2).*...
			 conv_results(:,:,:,:,3);
	T(:,:,:,:,2,1) = T(:,:,:,:,1,2);
	T(:,:,:,:,1,3) = (conv_results(:,:,:,:,6)+...
			  conv_results(:,:,:,:,8)).*...
			 conv_results(:,:,:,:,11)+...
			 conv_results(:,:,:,:,10).*...
			 conv_results(:,:,:,:,13)+...
			 conv_results(:,:,:,:,12).*...
			 conv_results(:,:,:,:,15)+...
			 gamma * conv_results(:,:,:,:,2).*...
			 conv_results(:,:,:,:,4);
	T(:,:,:,:,3,1) = T(:,:,:,:,1,3);
	T(:,:,:,:,1,4) = (conv_results(:,:,:,:,6)+...
			  conv_results(:,:,:,:,9)).*...
			 conv_results(:,:,:,:,12)+...
			 conv_results(:,:,:,:,10).*...
			 conv_results(:,:,:,:,14)+...
			 conv_results(:,:,:,:,11).*...
			 conv_results(:,:,:,:,15)+...
			 gamma * conv_results(:,:,:,:,2).*...
			 conv_results(:,:,:,:,5);
	T(:,:,:,:,4,1) = T(:,:,:,:,1,4);
	T(:,:,:,:,2,3) = (conv_results(:,:,:,:,7)+...
			  conv_results(:,:,:,:,8)).*...
			 conv_results(:,:,:,:,13)+...
			 conv_results(:,:,:,:,10).*...
			 conv_results(:,:,:,:,11)+...
			 conv_results(:,:,:,:,14).*...
			 conv_results(:,:,:,:,15)+...
			 gamma * conv_results(:,:,:,:,3).*...
			 conv_results(:,:,:,:,4);
	T(:,:,:,:,3,2) = T(:,:,:,:,2,3);
	T(:,:,:,:,2,4) = (conv_results(:,:,:,:,7)+...
			  conv_results(:,:,:,:,9)).*...
			 conv_results(:,:,:,:,14)+...
			 conv_results(:,:,:,:,10).*...
			 conv_results(:,:,:,:,12)+...
			 conv_results(:,:,:,:,13).*...
			 conv_results(:,:,:,:,15)+...
			 gamma * conv_results(:,:,:,:,3).*...
			 conv_results(:,:,:,:,5);
	T(:,:,:,:,4,2) = T(:,:,:,:,2,4);
	T(:,:,:,:,3,4) = (conv_results(:,:,:,:,8)+...
			  conv_results(:,:,:,:,9)).*...
			 conv_results(:,:,:,:,15)+...
			 conv_results(:,:,:,:,11).*...
			 conv_results(:,:,:,:,12)+...
			 conv_results(:,:,:,:,13).*...
			 conv_results(:,:,:,:,14)+...
			 gamma * conv_results(:,:,:,:,4).*...
			 conv_results(:,:,:,:,5);
	T(:,:,:,:,4,3) = T(:,:,:,:,3,4);
	
	T = squeeze(T);
    otherwise
	error('More than four dimensions are not supported.')
end

if nargout > 1
    params.spatial_size = spatial_size;
    params.region_of_interest = region_of_interest;
    params.gamma = gamma;
    params.sigma = sigma;
    params.delta = delta;
    params.c = certainty;
end
