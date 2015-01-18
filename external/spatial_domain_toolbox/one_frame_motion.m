function [d, c] = one_frame_motion(im, kernelsizes1, kernelsizes2, ...
				   model, method, d0)
% [D, C] = ESTIMATE_DISPLACEMENT(IM, KERNELSIZES1, KERNELSIZES2, ...
%                                MODEL, METHOD, D0)
%
% Estimate displacement from a single interlaced frame according to the
% algorithm described in section 8.6 of Gunnar Farnebäck's thesis
% "Polynomial Expansion for Orientation and Motion Estimation".
%
% IM
%   One interlaced image.
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

if nargin < 6
    d0 = zeros([size(im) 2]);
end

sides = size(im);

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
	    [A1, b1, c1, A2, b2, c2] = make_Abc_interlace(im, kernelsize1);

            % Ad hoc deweighting of expansion coefficients close to the border.
	    cin = ones(sides);
	    half = floor(kernelsize1/2);
	    border = get_border(sides, half);
	    cin(border) = 0.05 * cin(border);
	    border = get_border(sides, floor(half/2));
	    cin(border) = 0.05 * cin(border);
	else
	    % Ad hoc deweighting of expansion coefficients close to theborder.
	    cin = ones(sides);
	    half = floor(kernelsize1/2);
	    border = get_border(sides, half);
	    cin(border) = 0.2 * cin(border);
	    border = get_border(sides, floor(half/2));
	    cin(border) = 0.1 * cin(border);
	    
	    c_interlace = repmat(mod(1:sides(1), 2)', [1 sides(2)]);
	    
	    r1 = polyexp(im, c_interlace, 'quadratic', kernelsize1);
	    A1 = cat(4, cat(3, r1(:,:,4), r1(:,:,6)/2), ...
		     cat(3, r1(:,:,6)/2, r1(:,:,5)));
	    b1 = cat(3, r1(:,:,2), r1(:,:,3));
	    r2 = polyexp(im, (1-c_interlace), 'quadratic', kernelsize1);
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