function a = gaussian_app(side, dims, sigma)
% A = GAUSSIAN_APP(SIDE, DIMS, SIGMA)
%
% Construct a Gaussian applicability.
%
% SIDE    - Size of neighborhood along each dimension, should be odd.
% DIMS    - Number of signal dimensions, must be between 1 and 4.
% SIGMA   - Standard deviation of the Gaussian.
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

if nargin < 2
    error('missing argument')
end

if mod(side,2) == 0
    warning('side should be odd, increased it by one.')
    side = side + 1;
elseif mod(side,2) ~= 1
    error('side should be an odd integer.')
end

if side == 1
    a = 1;
    return
end

n = (side-1)/2;
I = (-n:n)';
switch dims
    case 1
	x = I;
	r = abs(x);
    case 2
	[x,y] = ndgrid(I);
	r = sqrt(x.^2+y.^2);
    case 3
	[x,y,z] = ndgrid(I);
	r = sqrt(x.^2+y.^2+z.^2);
    case 4
	[x,y,z,t] = ndgrid(I);
	r = sqrt(x.^2+y.^2+z.^2+t.^2);
end

if nargin < 3
    delta = 0.005;
    sigma = n/sqrt(-2*log(delta));
end

if sigma ~= 0
    a = exp(-r.^2/(2*sigma^2));
else
    a = double(r == 0);
end
