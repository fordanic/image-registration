function T = remove_isotropic(T)
% T = remove_isotropic(T)
%
% T is assumed to be an MxNx3x3 tensor field, where each 3x3 matrix is
% assumed to be symmetric. The function computes the smallest eigenvalue of
% each 3x3 matrix and subtracts this eigenvalue times the identity matrix.
% The algorithm follows appendix G of Gunnar Farnebäck's thesis "Polynomial
% Expansion for Orientation and Motion Estimation".
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

T(:,:,[1 5 9]) = T(:,:,[1 5 9]) - repmat(sum(T(:,:,[1 5 9]), 3)/3, [1 1 3]);
p = T(:,:,1,1).*T(:,:,2,2) + T(:,:,1,1).*T(:,:,3,3) + ...
    T(:,:,2,2).*T(:,:,3,3) - T(:,:,1,2).^2 - T(:,:,1,3).^2 - T(:,:,2,3).^2;
q = T(:,:,1,1).*T(:,:,2,3).^2 + T(:,:,2,2).*T(:,:,1,3).^2 + ...
    T(:,:,3,3).*T(:,:,1,2).^2 - 2*T(:,:,1,2).*T(:,:,1,3).*T(:,:,2,3) - ...
    T(:,:,1,1).*T(:,:,2,2).*T(:,:,3,3);
beta = sqrt(-4*p/3);

% eps is needed in the case that p=q=0. Notice that p<=0.
% Due to numerical errors it may happen that we get complex
% angles. It is correct to simply remove the imaginary part
% in these cases.
alpha = real(acos(3*q./(-eps+p.*beta))/3);
lambda3 = beta.*cos(alpha+2*pi/3);
T(:,:,[1 5 9]) = T(:,:,[1 5 9]) - repmat(lambda3, [1 1 3]);
