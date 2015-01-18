function p = min_quadform(Qtot)
% P = MIN_QUADFORM(Qtot)
% 
% Minimize quadratic forms according to equations 6.20 -- 6.24 in
% Gunnar Farnebäck's thesis "Polynomial Expansion for Orientation and
% Motion Estimation".
%
% Qtot   - A collection of quadratic forms, having the size
%          HEIGTH x WIDTH x N x N.
% 
% P      - A collection of optimal parameters, having the size
%          HEIGHT x WIDTH x (N-1).
% 
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

error('MIN_QUADFORM is implemented as a mex-file. It has not been compiled on this platform.')
