function [d, c] = estimate_displacement(im1, im2, kernelsizes1, ...
					kernelsizes2, model, method, d0)
% [D, C] = ESTIMATE_DISPLACEMENT(IM1, IM2, KERNELSIZES1, KERNELSIZES2, ...
%                                MODEL, METHOD, D0)
%
% Estimate displacement according to the algorithms described in chapter 7
% of Gunnar Farnebäck's thesis "Polynomial Expansion for Orientation and
% Motion Estimation".
%
% IM1, IM2
%   Two grayscale images of the same size.
% 
% KERNELSIZES1
%   Vector of kernelsizes used for polynomial expansion in each iteration
%   of the displacement computations. The same kernelsize can be repeated.
%
% KERNELSIZES2    
%   Vector of kernelsizes used for averaging each iteration of the
%   displacement computations. Must have the same length as KERNELSIZES1.
% 
% MODEL
%   Choice of parametric motion model, 'constant', 'affine', or 'eightparam'.
%
% METHOD
%   'fast' or 'accurate'
%
% D0 [optional]
%   A priori displacement estimate. Default is an all-zero displacement
%   field.
%
% D
%   Estimated displacement field.
%
% C
%   Estimated (reversed) confidence value. Small values indicate more
%   reliable displacement estimates.
%
%
% Warning: This function may be revised in a non-backwards compatible
%          way in later revisions of the Spatial Domain Toolbox.
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se
   
if length(kernelsizes1) ~= length(kernelsizes2)
    error('kernelsizes1 and kernelsizes2 must have the same length')
end

if nargin < 7
    d0 = zeros([size(im1) 2]);
end

d = d0;

for k = 1:length(kernelsizes1)
    kernelsize1 = kernelsizes1(k);
    kernelsize2 = kernelsizes2(k);
    
    if k == 1
	last_kernelsize1 = -1;
    else
	last_kernelsize1 = kernelsizes1(k - 1);
    end
    
    if kernelsize1 ~= last_kernelsize1
	if strcmp(method, 'fast')
	    [A1, b1, c1] = make_Abc_fast(im1, kernelsize1);
	    [A2, b2, c2] = make_Abc_fast(im2, kernelsize1);

            % Ad hoc deweighting of expansion coefficients close to the border.
	    cin = ones(size(im1));
	    half = floor(kernelsize1/2);
	    border = get_border(size(cin), half);
	    cin(border) = 0.05 * cin(border);
	    border = get_border(size(cin), floor(half/2));
	    cin(border) = 0.05 * cin(border);
	else
	    % Ad hoc deweighting of expansion coefficients close to theborder.
	    cin = ones(size(im1));
	    half = floor(kernelsize1/2);
	    border = get_border(size(cin), half);
	    cin(border) = 0.2 * cin(border);
	    border = get_border(size(cin), floor(half/2));
	    cin(border) = 0.1 * cin(border);
	    
	    r1 = polyexp(im1, cin, 'quadratic', kernelsize1);
	    A1 = cat(4, cat(3, r1(:,:,4), r1(:,:,6)/2), ...
		     cat(3, r1(:,:,6)/2, r1(:,:,5)));
	    b1 = cat(3, r1(:,:,2), r1(:,:,3));
	    r2 = polyexp(im2, cin, 'quadratic', kernelsize1);
	    A2 = cat(4, cat(3, r2(:,:,4), r2(:,:,6)/2), ...
		     cat(3, r2(:,:,6)/2, r2(:,:,5)));
	    b2 = cat(3, r2(:,:,2), r2(:,:,3));
	end
    end
    
    d0 = d;
    sigma = 0.15 * (kernelsize2 - 1);
    [A, b] = prepare_displacement_matrices(A1, b1, A2, b2, d0);
    [d, c] = compute_displacement(A, b, kernelsize2, sigma, cin, model);
end

return

function border = get_border(sides, width)

border = zeros(sides);
border(1:width, :) = 1;
border(end-width+1:end, :) = 1;
border(:, 1:width) = 1;
border(:, end-width+1:end) = 1;
border = logical(border);

return