function d = estimate_disparity(leftim, rightim, kernelsize1, sigma1, ...
				kernelsize2, sigma2, dmax)
    
% D = ESTIMATE_DISPARITY(LEFTIM, RIGHTIM, KERNELSIZE1, SIGMA1, ...
%                        KERNELSIZE2, SIGMA2, DMAX)
%
% Estimate disparity from a stereo pair, using the algorithm described
% in section 7.3 of Gunnar Farnebäck's thesis "Polynomial Expansion for
% Orientation and Motion Estimation".
%
% LEFTIM
%   Left image of stereo pair.
%
% RIGHTIM
%   Left image of stereo pair.
%
% KERNELSIZE1
%   Size of Gaussian applicability used for polynomial expansion.
%
% SIGMA1
%   Standard deviation of Gaussian applicability used for polynomial
%   expansion.
%
% KERNELSIZE2
%   Size of Gaussian applicability used for normalized averaging.
%
% SIGMA2
%   Standard deviation of Gaussian applicability used for normalized
%   averaging.
%
% DMAX
%   Upper bound on the disparity values that can occur in the stereo pair.
% 
% D
%   Estimated disparity values.
%
%
% Note: It is recommended to have kernelsize2 > kernelsize1.
%       Choose kernelsize1 as small as possible without getting
%       a significantly anisotropic applicability. 
%
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se
    
options.sigma = sigma1;
[Al, bl, cl] = make_Abc_fast(leftim, kernelsize1, [], options);
[Ar, br, cr] = make_Abc_fast(rightim, kernelsize1, [], options);

A = (Al + Ar) / 2;
b = -(bl - br) / 2;

% Vectorized computation of equation (7.17).
det = A(:,:,1,1).*A(:,:,2,2) - A(:,:,1,2).^2;
dy = ( A(:,:,2,2).*b(:,:,1) - A(:,:,1,2).*b(:,:,2)) ./ det;
dx = (-A(:,:,2,1).*b(:,:,1) + A(:,:,1,1).*b(:,:,2)) ./ det;

c1 = dx.^2 ./ (eps + dy.^2 + dx.^2);

% Outlier rejection.
c2 = ones(size(leftim));
c2(dx <  0) = 0;
c2(dx > dmax) = 0;

% Set certainty to zero along the borders. The polynomial expansion is not
% reliable there.
M = (kernelsize1 - 1) / 2;
c3 = ones(size(leftim));
c3(1:M,:)         = 0;
c3(end-M+1:end,:) = 0;
c3(:,1:M)         = 0;
c3(:,end-M+1:end) = 0;

certainty = c1 .* c2 .* c3;
x = -(kernelsize2-1)/2 : (kernelsize2-1)/2;
a = exp(-x.^2/(2*sigma2^2));

% Normalized averaging, computed with separable correlations.
d = conv3(conv3(certainty.*dx, a), a') ./ conv3(conv3(certainty, a), a'); 
