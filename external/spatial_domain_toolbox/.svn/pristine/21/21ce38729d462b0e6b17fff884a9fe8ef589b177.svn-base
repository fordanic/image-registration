function [displacement, cout] = compute_displacement(A, b, kernelsize, ...
						     sigma, cin, model)
% [DISPLACEMENT, COUT] = COMPUTE_DISPLACEMENT(A, B, KERNELSIZE, SIGMA, ...
%                                             CIN, MODEL)
%
% Compute displacement estimates according to equation (7.30) in Gunnar
% Farnebäck's thesis "Polynomial Expansion for Orientation and Motion
% Estimation". Optionally also compute output certainty according to
% equation (7.24) extended to parameterized displacement fields.
%
% A, B
%   Displacement matrices computed by PREPARE_DISPLACEMENT_MATRICES
%
% KERNELSIZE, SIGMA
%   Size and standard deviation for the Gaussian applicability used in
%   averaging.
%
% CIN
%   Input certainty.
%
% MODEL
%   Choice of parametric motion model, 'constant', 'affine', or 'eightparam'.
%
% DISPLACEMENT
%   Computed displacement field.
%
% COUT
%   Computed (reversed) confidence value. Small values indicate more
%   reliable displacement values.
%
% Warning: This function may be revised in a non-backwards compatible
%          way in later revisions of the Spatial Domain Toolbox.
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

sides = size(A);
sides = sides(1:2);

switch model
 case 'constant'

  Q = zeros([sides 5]);
  Q(:,:,1) = A(:,:,1,1).^2 + A(:,:,1,2).^2;
  Q(:,:,2) = A(:,:,2,2).^2 + A(:,:,1,2).^2;
  Q(:,:,3) = (A(:,:,1,1) + A(:,:,2,2)).*A(:,:,1,2);
  Q(:,:,4) = A(:,:,1,1).*b(:,:,1) + A(:,:,1,2).*b(:,:,2);
  Q(:,:,5) = A(:,:,1,2).*b(:,:,1) + A(:,:,2,2).*b(:,:,2);

  app = gaussian_app(kernelsize, 1, sigma);
  cinaver = conv3(conv3(cin, app), app');
  Q = conv3(conv3(Q.*repmat(cin, [1 1 5]), app), app') ./ ...
      (eps + repmat(cinaver, [1 1 5])); 
  
  % Solve the equation Qv=q.
  a = Q(:,:,1);
  bb = Q(:,:,2);
  c = Q(:,:,3);
  d = Q(:,:,4);
  e = Q(:,:,5);
  % Q=[a c;c b], inv(Q)=[b -c;-c a]/(a*b-c^2), q=[d;e]
  displacement = cat(3,d.*bb-c.*e,a.*e-c.*d)./(eps+repmat(a.*bb-c.^2,[1 1 2]));
  
  % Compute output certainty
  if nargout > 1
      q = b(:,:,1).^2 + b(:,:,2).^2;
      q = conv3(conv3(q.*cin, app), app') ./ (eps + cinaver); 
      cout = q - d.*displacement(:,:,1) - e.*displacement(:,:,2);
  end
  
 case 'affine'
  [x,y] = ndgrid(1:sides(1), 1:sides(2));
%  [x,y] = meshgrid(1:sides(2), 1:sides(1));
  Q = zeros([sides 24]);
  Q(:,:,1)  = A(:,:,1,1).^2 + A(:,:,1,2).^2;                % (1,1)
  Q(:,:,2)  = Q(:,:,1).*x;                                  % (1,2) (2,1)
  Q(:,:,3)  = Q(:,:,1).*y;                                  % (1,3) (3,1)
  Q(:,:,4)  = (A(:,:,1,1) + A(:,:,2,2)).*A(:,:,1,2);        % (1,4) (4,1)
  Q(:,:,5)  = Q(:,:,4).*x;                      % (1,5) (5,1) (2,4) (4,2)
  Q(:,:,6)  = Q(:,:,4).*y;                      % (1,6) (6,1) (3,4) (4,3)
  Q(:,:,7)  = Q(:,:,2).*x;				    % (2,2)
  Q(:,:,8)  = Q(:,:,2).*y;				    % (2,3) (3,2)
  Q(:,:,9)  = Q(:,:,5).*x;				    % (2,5) (5,2)
  Q(:,:,10) = Q(:,:,5).*y;                      % (2,6) (6,2) (3,5) (5,3)
  Q(:,:,11) = Q(:,:,3).*y;				    % (3,3)
  Q(:,:,12) = Q(:,:,6).*y;				    % (3,6) (6,3)
  Q(:,:,13) = A(:,:,1,2).^2 + A(:,:,2,2).^2;                % (4,4)
  Q(:,:,14) = Q(:,:,13).*x;				    % (4,5) (5,4)
  Q(:,:,15) = Q(:,:,13).*y;				    % (4,6) (6,4)
  Q(:,:,16) = Q(:,:,14).*x;				    % (5,5)
  Q(:,:,17) = Q(:,:,14).*y;				    % (5,6) (6,5)
  Q(:,:,18) = Q(:,:,15).*y;				    % (6,6)

  Q(:,:,19) = A(:,:,1,1).*b(:,:,1) + A(:,:,1,2).*b(:,:,2);  % (1)
  Q(:,:,20) = Q(:,:,19).*x;				    % (2)
  Q(:,:,21) = Q(:,:,19).*y;				    % (3)
  Q(:,:,22) = A(:,:,1,2).*b(:,:,1) + A(:,:,2,2).*b(:,:,2);  % (4)
  Q(:,:,23) = Q(:,:,22).*x;				    % (5)
  Q(:,:,24) = Q(:,:,22).*y;				    % (6)


  % Compute displacement from affine fields in each neighborhood.
  app = gaussian_app(kernelsize, 1, sigma);
  cinaver = conv3(conv3(cin, app), app');
  Q = conv3(conv3(Q.*repmat(cin, [1 1 24]), app), app') ./ ...
      (eps + repmat(cinaver, [1 1 24]));
  % We build the equation system as a quadratic form that min_quadform
  % can solve. Slightly wasteful but effective.
  Q = reshape(Q(:,:,...
		[ 1  2  3  4  5  6 19
		  2  7  8  5  9 10 20
		  3  8 11  6 10 12 21
		  4  5  6 13 14 15 22
		  5  9 10 14 16 17 23
		  6 10 12 15 17 18 24
		 19 20 21 22 23 24 24]),[sides 7 7]);

  % Solve the equation Qv=q.
  p = -min_quadform(Q);
  displacement = zeros([sides 2]);
  displacement(:,:,1) = sum(p(:,:,1:3).*cat(3, ones(sides), x, y), 3);
  displacement(:,:,2) = sum(p(:,:,4:6).*cat(3, ones(sides), x, y), 3);

  if nargout > 1
      q = b(:,:,1).*b(:,:,1) + b(:,:,2).*b(:,:,2);
      q = conv3(conv3(q.*cin, app), app') ./ (eps + cinaver);
      cout = q - sum(p .* Q(:,:,1:6,7), 3);
  end
  
 case 'eightparam'
  [x,y] = ndgrid(1:sides(1), 1:sides(2));
  Q = zeros([sides 39]);
  Q(:,:,1)  = A(:,:,1,1).^2 + A(:,:,1,2).^2;                % (1,1)
  Q(:,:,2)  = Q(:,:,1).*x;                                  % (1,2) (2,1)
  Q(:,:,3)  = Q(:,:,1).*y;                                  % (1,3) (3,1)
  Q(:,:,4)  = (A(:,:,1,1) + A(:,:,2,2)).*A(:,:,1,2);        % (1,4) (4,1)
  Q(:,:,5)  = Q(:,:,4).*x;                      % (1,5) (5,1) (2,4) (4,2)
  Q(:,:,6)  = Q(:,:,4).*y;                      % (1,6) (6,1) (3,4) (4,3)
  Q(:,:,7)  = Q(:,:,2).*x;				    % (2,2)
  Q(:,:,8)  = Q(:,:,2).*y;				    % (2,3) (3,2)
  Q(:,:,9)  = Q(:,:,5).*x;				    % (2,5) (5,2)
  Q(:,:,10) = Q(:,:,5).*y;                      % (2,6) (6,2) (3,5) (5,3)
  Q(:,:,11) = Q(:,:,3).*y;				    % (3,3)
  Q(:,:,12) = Q(:,:,6).*y;				    % (3,6) (6,3)
  Q(:,:,13) = A(:,:,1,2).^2 + A(:,:,2,2).^2;                % (4,4)
  Q(:,:,14) = Q(:,:,13).*x;				    % (4,5) (5,4)
  Q(:,:,15) = Q(:,:,13).*y;				    % (4,6) (6,4)
  Q(:,:,16) = Q(:,:,14).*x;				    % (5,5)
  Q(:,:,17) = Q(:,:,14).*y;				    % (5,6) (6,5)
  Q(:,:,18) = Q(:,:,15).*y;				    % (6,6)
  Q(:,:,19) = Q(:,:,7) + Q(:,:,10);			    % (1,7) (7,1)
  Q(:,:,20) = Q(:,:,19).*x;				    % (2,7) (7,2)
  Q(:,:,21) = Q(:,:,19).*y;		        % (3,7) (7,3) (2,8) (8,2)
  Q(:,:,22) = Q(:,:,9) + Q(:,:,17);			    % (4,7) (7,4)
  Q(:,:,23) = Q(:,:,22).*x;				    % (5,7) (7,5)
  Q(:,:,24) = Q(:,:,22).*y;		        % (6,7) (7,6) (5,8) (8,5)
  Q(:,:,25) = Q(:,:,8) + Q(:,:,12);			    % (1,8) (8,1)
  Q(:,:,26) = Q(:,:,25).*y;				    % (3,8) (8,3)
  Q(:,:,27) = Q(:,:,10) + Q(:,:,18);			    % (4,8) (8,4)
  Q(:,:,28) = Q(:,:,27).*y;				    % (6,8) (8,6)
  Q(:,:,29) = (Q(:,:,20) + Q(:,:,24)).*x;		    % (7,7)
  Q(:,:,30) = (Q(:,:,21) + Q(:,:,28)).*x;		    % (7,8) (8,7)
  Q(:,:,31) = (Q(:,:,21) + Q(:,:,28)).*y;		    % (8,8)
  							     
  Q(:,:,32) = A(:,:,1,1).*b(:,:,1) + A(:,:,1,2).*b(:,:,2);  % (1)
  Q(:,:,33) = Q(:,:,32).*x;				    % (2)
  Q(:,:,34) = Q(:,:,32).*y;				    % (3)
  Q(:,:,35) = A(:,:,1,2).*b(:,:,1) + A(:,:,2,2).*b(:,:,2);  % (4)
  Q(:,:,36) = Q(:,:,35).*x;				    % (5)
  Q(:,:,37) = Q(:,:,35).*y;				    % (6)
  Q(:,:,38) = (Q(:,:,33) + Q(:,:,37)).*x;		    % (7)
  Q(:,:,39) = (Q(:,:,33) + Q(:,:,37)).*y;                   % (8)


  % Compute displacement from eightparam fields in each neighborhood.
  app = gaussian_app(kernelsize, 1, sigma);
  cinaver = conv3(conv3(cin, app), app');
  Q = conv3(conv3(Q.*repmat(cin, [1 1 39]), app), app') ./ ...
      (eps + repmat(cinaver, [1 1 39]));
  % We build the equation system as a quadratic form that min_quadform
  % can solve. Slightly wasteful but effective.
  Q = reshape(Q(:,:,...
		[ 1  2  3  4  5  6 19 25 32
		  2  7  8  5  9 10 20 21 33
		  3  8 11  6 10 12 21 26 34
		  4  5  6 13 14 15 22 27 35
		  5  9 10 14 16 17 23 24 36
		  6 10 12 15 17 18 24 28 37
		 19 20 21 22 23 24 29 30 38
		 25 21 26 27 24 28 30 31 39
		 32 33 34 35 36 37 38 39 39]),[sides 9 9]);

  % Solve the equation Qv=-q.
  p = -min_quadform(Q);
  displacement = zeros([sides 2]);
  displacement(:,:,1) = sum(p(:,:,[1 2 3 7 8]).*cat(3, ones(sides), ...
						    x, y, x.^2, x.*y), 3);
  displacement(:,:,2) = sum(p(:,:,4:8).*cat(3, ones(sides), ...
					    x, y, x.*y, y.^2), 3);

  if nargout > 1
      q = b(:,:,1).*b(:,:,1) + b(:,:,2).*b(:,:,2);
      q = conv3(conv3(q.*cin, app), app') ./ (eps + cinaver);
      cout = q - sum(p .* Q(:,:,1:8,9), 3);
  end
end
