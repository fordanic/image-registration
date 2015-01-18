function [b, c, params] = make_bc_fast(signal, spatial_size, options)
% MAKE_BC_FAST Performs a linear polynomial expansion
%
% [b, c] = make_bc_fast(signal, spatial_size, options)
% or
% [b, c, params] = make_bc_fast(signal, spatial_size, options)
%
% Compute b, and c parameters in up to four dimensions. The
% parameters relate to the local signal model
% f(x) = b^T x + c
% and are determined by a Gaussian weighted least squares fit. This
% implementation uses a fast hierarchical scheme of separable filters,
% described in chapter 4 of "Polynomial Expansion for Orientation and
% Motion Estimation" by Gunnar Farnebäck.
%
% signal                        - Signal values. Must be real and nonsparse
%                                 and the number of dimensions, N, must be
%                                 at most four.
%
% spatial_size [optional]       - Size of the spatial support of the filters
%                                 along each dimension. Default value is 9.
%
% options [optional]            - Struct array that may contain various
%                                 parameters that affect the algorithm.
%                                 These are explained below.
%
% A                             - Computed A matrices. A has N+2
%                                 dimensions, where the first N indices
%                                 indicates the position in the signal and
%                                 the last two contains the matrix for each
%                                 point. In the case that region_of_interest
%                                 is less than N-dimensional, the singleton
%                                 dimensions are removed.
%
% b                             - Computed b vectors. b has N+1
%                                 dimensions, where the first N indices
%                                 indicates the position in the signal and
%                                 the last one contains the vector for each
%                                 point. In the case that region_of_interest
%                                 is less than N-dimensional, the singleton
%                                 dimensions are removed.
%
% c                             - Computed c scalars. c has N dimensions.
%                                 In the case that region_of_interest
%                                 is less than N-dimensional, the singleton
%                                 dimensions are removed.
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
% This is an adapted version of make_Abc_fast
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

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

N = ndims(signal);
if N == 2 && size(signal, 2) == 1
    N = 1;
end

if nargin < 2
    spatial_size = 9;
end

if spatial_size < 1
    error('What use would such a small kernel be?')
elseif mod(spatial_size, 2) ~= 1
    spatial_size = 2*floor((spatial_size-1)/2)+1;
    warning(sprintf('Only kernels of odd size are allowed. Changed the size to %d.', spatial_size))
end

sigma = 0.15 * (spatial_size - 1);
certainty = ones([repmat(spatial_size, [1 N]) 1]);
n = (spatial_size-1)/2;

if nargin == 3
    if isfield(options, 'sigma')
        sigma = options.sigma;
    elseif isfield(options, 'delta')
        sigma = n/sqrt(-2*log(delta));
    end
    if isfield(options, 'c')
        certainty = options.c;
    end
end

a = exp(-(-n:n).^2/(2*sigma^2))';

switch N
    case 1
        % Set up applicability and basis functions.
        applicability = a;
        x = (-n:n)';
        b = [ones(size(x)) x];
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
        conv_results = zeros([size(signal) 2]);
        conv_results(:,1) = imfilter(signal, kernelx0);
        conv_results(:,2) = imfilter(signal, kernelx1);
        
        % Apply the inverse metric.
        conv_results(:,1) = Qinv(1,1)*conv_results(:,1);
        conv_results(:,2) = Qinv(2,2)*conv_results(:,2);
        
        % Build b, and c.
        %
        % b, and c are obtained from the convolution results according to
        %
        % b=[2], c=[1].
        %
        
        b = zeros([size(signal) 1]);
        b(:,1) = conv_results(:,2);
        c = conv_results(:,1);
        
        b = squeeze(b);
        c = squeeze(c);
        
    case 2
        % Set up applicability and basis functions.
        applicability = a*a';
        [x,y] = ndgrid(-n:n);
        b = cat(3, ones(size(x)), x, y);
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
        conv_y0 = imfilter(signal, kernely0);
        conv_y1 = imfilter(signal, kernely1);
        
        % Convolutions in the x-direction.
        kernelx0 = kernely0(:);
        kernelx1 = kernely1(:);
        conv_results = zeros([size(signal) 3]);
        conv_results(:,:,1) = imfilter(conv_y0, kernelx0); % y0x0
        conv_results(:,:,2) = imfilter(conv_y0, kernelx1); % y0x1
        clear conv_y0
        conv_results(:,:,3) = imfilter(conv_y1, kernelx0); % y1x0
        clear conv_y1
        
        % Apply the inverse metric.
        conv_results(:,:,1) = Qinv(1,1)*conv_results(:,:,1);
        conv_results(:,:,2) = Qinv(2,2)*conv_results(:,:,2);
        conv_results(:,:,3) = Qinv(3,3)*conv_results(:,:,3);
        
        % Build b, and c.
        %
        % b, and c are obtained from the convolution results according to
        %
        %   [2]
        % b=[3], c=[1].
        %
        % where the off-diagonal elements in A additionally are halved.
        %
        
        b = zeros([size(signal) 2]);
        b(:,:,1) = conv_results(:,:,2);
        b(:,:,2) = conv_results(:,:,3);
        c = conv_results(:,:,1);
        
        b = squeeze(b);
        c = squeeze(c);
        
    case 3
        % Set up applicability and basis functions.
        applicability = outerprod(a, a, a);
        [x,y,t] = ndgrid(-n:n);
        b = cat(4, ones(size(x)), x, y, t);
        nb = size(b,4);
        
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
        kernelt0 = reshape(a, [1 1 spatial_size]);
        kernelt1 = reshape((-n:n)'.*a, [1 1 spatial_size]);
        conv_t0 = imfilter(signal, kernelt0);
        conv_t1 = imfilter(signal, kernelt1);
        
        % Convolutions in the y-direction
        kernely0 = reshape(kernelt0, [1 spatial_size]);
        kernely1 = reshape(kernelt1, [1 spatial_size]);
        conv_t0y0 = imfilter(conv_t0, kernely0);
        conv_t0y1 = imfilter(conv_t0, kernely1);
        clear conv_t0
        conv_t1y0 = imfilter(conv_t1, kernely0);
        
        % Convolutions in the x-direction
        kernelx0 = reshape(kernelt0, [spatial_size 1]);
        kernelx1 = reshape(kernelt1, [spatial_size 1]);
        conv_results = zeros([size(signal) 4]);
        conv_results(:,:,:,1) = imfilter(conv_t0y0, kernelx0); % t0y0x0
        conv_results(:,:,:,2) = imfilter(conv_t0y0, kernelx1); % t0y0x1
        clear conv_t0y0
        conv_results(:,:,:,3) = imfilter(conv_t0y1, kernelx0); % t0y1x0
        clear conv_t0y1
        conv_results(:,:,:,4) = imfilter(conv_t1y0, kernelx0); % t1y0x0
        clear conv_t1y0
        
        % Apply the inverse metric.
        conv_results(:,:,:,1) = Qinv(1,1)*conv_results(:,:,:,1);
        conv_results(:,:,:,2) = Qinv(2,2)*conv_results(:,:,:,2);
        conv_results(:,:,:,3) = Qinv(3,3)*conv_results(:,:,:,3);
        conv_results(:,:,:,4) = Qinv(4,4)*conv_results(:,:,:,4);
        
        % Build b, and c.
        %
        % b, and c are obtained from the convolution results according to
        %
        %   [2]
        % b=[3], c=[1]
        %   [4]
        %
        % where the off-diagonal elements in A additionally are halved.
        %
        
        b = zeros([size(signal) 3]);
        b(:,:,:,1) = conv_results(:,:,:,2);
        b(:,:,:,2) = conv_results(:,:,:,3);
        b(:,:,:,3) = conv_results(:,:,:,4);
        c = conv_results(:,:,:,1);
        
        b = squeeze(b);
        c = squeeze(c);
        
    case 4
        % Set up applicability and basis functions.
        applicability = outerprod(a, a, a, a);
        [x,y,z,t] = ndgrid(-n:n);
        b = cat(5, ones(size(x)), x, y, z, t);
        nb = size(b, 5);
        
        % Compute the inverse metric.
        Q = zeros(nb, nb);
        for i = 1:nb
            for j = i:nb
                Q(i,j) = sum(sum(sum(sum(b(:,:,:,:,i).*applicability.*certainty.*b(:,:,:,:,j)))));
                Q(j,i) = Q(i,j);
            end
        end
        clear b applicability x y z t
        Qinv = inv(Q);
        
        % Convolutions in the t-direction
        kernelt0 = reshape(a, [1 1 1 spatial_size]);
        kernelt1 = reshape((-n:n)'.*a, [1 1 1 spatial_size]);
        conv_t0 = imfilter(signal, kernelt0);
        conv_t1 = imfilter(signal, kernelt1);
        
        % Convolutions in the z-direction
        kernelz0 = reshape(kernelt0, [1 1 spatial_size]);
        kernelz1 = reshape(kernelt1, [1 1 spatial_size]);
        conv_t0z0 = imfilter(conv_t0, kernelz0);
        conv_t0z1 = imfilter(conv_t0, kernelz1);
        clear conv_t0
        conv_t1z0 = imfilter(conv_t1, kernelz0);
        clear conv_t1
        
        % Convolutions in the y-direction
        kernely0 = reshape(kernelt0, [1 spatial_size]);
        kernely1 = reshape(kernelt1, [1 spatial_size]);
        conv_t0z0y0 = imfilter(conv_t0z0, kernely0);
        conv_t0z0y1 = imfilter(conv_t0z0, kernely1);
        clear conv_t0z0
        conv_t0z1y0 = imfilter(conv_t0z1, kernely0);
        clear conv_t0z1
        conv_t1z0y0 = imfilter(conv_t1z0, kernely0);
        clear conv_t1z0
        
        % Convolutions in the x-direction
        kernelx0 = reshape(kernelt0, [spatial_size 1]);
        kernelx1 = reshape(kernelt1, [spatial_size 1]);
        conv_results = zeros([size(signal) 5]);
        conv_results(:,:,:,:,1) = imfilter(conv_t0z0y0, kernelx0); % t0z0y0x0
        conv_results(:,:,:,:,2) = imfilter(conv_t0z0y0, kernelx1); % t0z0y0x1
        clear conv_t0z0y0
        conv_results(:,:,:,:,3) = imfilter(conv_t0z0y1, kernelx0); % t0z0y1x0
        clear conv_t0z0y1
        conv_results(:,:,:,:,4) = imfilter(conv_t0z1y0, kernelx0); % t0z1y0x0
        clear conv_t0z1y0
        conv_results(:,:,:,:,5) = imfilter(conv_t1z0y0, kernelx0); % t1z0y0x0
        clear conv_t1z0y0
        
        % Apply the inverse metric.
        conv_results(:,:,:,:,1) = Qinv(1,1)*conv_results(:,:,:,:,1);
        conv_results(:,:,:,:,2) = Qinv(2,2)*conv_results(:,:,:,:,2);
        conv_results(:,:,:,:,3) = Qinv(3,3)*conv_results(:,:,:,:,3);
        conv_results(:,:,:,:,4) = Qinv(4,4)*conv_results(:,:,:,:,4);
        conv_results(:,:,:,:,5) = Qinv(5,5)*conv_results(:,:,:,:,5);
        
        % Build b, and c.
        %
        % b, and c are obtained from the convolution results according to
        %
        %
        %   [2]
        %   [3]
        % b=[4]. c=[1]
        %   [5]
        % where the off-diagonal elements in A additionally are halved.
        %
        
        b = zeros([size(signal) 4]);
        b(:,:,:,:,1) = conv_results(:,:,:,:,2);
        b(:,:,:,:,2) = conv_results(:,:,:,:,3);
        b(:,:,:,:,3) = conv_results(:,:,:,:,4);
        b(:,:,:,:,4) = conv_results(:,:,:,:,5);
        c = conv_results(:,:,:,:,1);
        
        b = squeeze(b);
        c = squeeze(c);
        
        
    otherwise
        error('More than four dimensions are not supported.')
end

if nargout > 3
    params.spatial_size = spatial_size;
    params.region_of_interest = region_of_interest;
    params.sigma = sigma;
    params.delta = delta;
    params.c = certainty;
end
