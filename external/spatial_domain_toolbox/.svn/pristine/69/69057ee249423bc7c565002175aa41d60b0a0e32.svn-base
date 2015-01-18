function phi = angular_error(u, v)
% PHI = ANGULAR_ERROR(U, V) 
%
% Compute the pointwise angular error between two velocity fields. U and V
% are supposed to be arrays of dimensions MxNx2. PHI is given in degrees.
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

m = size(u, 1);
n = size(u, 2);

u3d = cat(3, u, ones(m,n));
v3d = cat(3, v, ones(m,n));

u3d_normalized = u3d./repmat(sqrt(sum(u3d.^2, 3)), [1 1 3]);
v3d_normalized = v3d./repmat(sqrt(sum(v3d.^2, 3)), [1 1 3]);

cosphi = sum(u3d_normalized.*v3d_normalized, 3);
phi = 180/pi*acos(cosphi);
