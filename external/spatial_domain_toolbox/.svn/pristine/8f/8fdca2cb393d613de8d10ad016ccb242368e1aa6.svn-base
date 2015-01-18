function [A, b] = prepare_displacement_matrices(A1, b1, A2, b2, displacement)

% [A, b] = PREPARE_DISPLACEMENT_MATRICES(A1, b1, A2, b2, displacement)
%
% Compute matrices used for displacement estimation as defined by equations
% (7.32) and (7.33) in Gunnar Farnebäck's thesis "Polynomial Expansion for
% Orientation and Motion Estimation".
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

% The code below has been replaced by a mex function with exactly the same
% functionality. Since this implementation is so much slower we give an
% error if the mex-file is missing.
    
error('PREPARE_DISPLACEMENT_MATRICES is implemented as a mex-file. It has not been compiled on this platform.')
    
sides = size(A1);
sides = sides(1:2);

if nargin < 5
    displacement = zeros([sides 2]);
end

A = zeros(size(A1));
b = zeros(size(b1));

% If displacement is zero, we will get A = (A1+A2)/2 and b = -(b2-b1)/2.

for j = 1:sides(2)
    for i = 1:sides(1)
	di = displacement(i,j,1);
	if i + di < 1
	    di = 1-i;
	end
	if i + di > sides(1)
	    di = sides(1) - i;
	end
	dj = displacement(i,j,2);
	if j + dj < 1
	    dj = 1-j;
	end
	if j + dj > sides(2)
	    dj = sides(2) - j;
	end
	
	A(i,j,:,:) = (A1(i,j,:,:) + A2(i+di,j+dj,:,:)) / 2;
	AA = squeeze(A(i,j,:,:));
	bb2 = squeeze(b2(i+di,j+dj,:)) - 2 * AA * [di;dj];
	b(i,j,:) = -(shiftdim(bb2,-2) - b1(i,j,:)) / 2;
    end
end
