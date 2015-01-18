function [phi, phi2] = evaluate_velocity(v, v0, mask, c, densities)
% [PHI, PHI2] = EVALUATE_VELOCITY(V, V0, MASK, C, DENSITIES)
%
% Evaluate velocity estimates, using the average spatiotemporal angular
% error.
%
% V    - Estimated velocity field
% V0   - Correct velocity field
% MASK - Ignore values where mask is zero
% C    - Confidence values. Partial coverage errors are computed using
%        the desired fraction of the pixels with smallest (!) confidence
%        values.
% DENSITIES - Set of density values to report. Default is [1.0, 0.9, 0.7].
%
% PHI  - Angular errors at each pixel, measured in degrees.
% PHI2 - Angular errors at each pixel under MASK, drawn out to a vector.
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

if ~islogical(mask)
    error('mask must be a logical array');
end

if nargin < 5
    densities = [1 0.9 0.7];
end

if nargin < 4
    c = ones(size(v));
    densities = 1;
end

phi = angular_error(v0, v);

phi2 = phi(mask);
c2 = c(mask);
[tmp, order] = sort(c2);
disp(sprintf('average   std    density'))
disp(sprintf(' error    dev'))
for density = densities
    m = mean(phi2(order(1:floor(density*length(phi2)))));
    s = std(phi2(order(1:floor(density*length(phi2)))));
    disp(sprintf('%6.3f   %6.3f  %3d%%', m, s, round(density*100)))
end
