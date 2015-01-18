function T = poly_to_tensor(r, gamma)
% T = POLY_TO_TENSOR(R, GAMMA)
%
% Construct orientation tensors in up to four dimensions from quadratic
% polynomial expansion coefficients. The algorithm is described in chapter
% 5 of Gunnar Farnebäck's thesis "Polynomial Expansion for Orientation and
% Motion Estimation".
%
% R       - Quadratic polynomial expansion coefficients, ordered as
%           in the result of the POLYEXP function. The last dimension
%           must be 3 for 1D, 6 for 2D, 10 for 3D, and 15 for 4D tensors.
%           The number of preceding dimensions in R does not matter.
%
% GAMMA   - Relation between the contribution to the tensor from the
%           linear and the quadratic parts of the signal, as specified
%           in equation (5.19). 0 means that only the quadratic part
%           matters while a very large number means that only the
%           linear part is used. A heuristically recommended value is
%           1/(8*sigma^2), where sigma is the standard deviation of the
%           applicability used in the polynomial expansion. Gamma can
%           also be specified as either of the strings 'even' and 'odd'.
%           The former is equivalent to setting gamma=0 while the latter
%           is the limit value of setting gamma very large and scaling
%           by 1/gamma. Using the value inf is equivalent to 'odd'.
%
% T       - Computed tensors. These have one dimension more than R and
%           the last dimensions are NxN, where N is the dimensionality of
%           the tensors. The size of the preceding dimensions are the
%           same as for R.
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

sides = size(r);
basis_size = sides(end);
sides = sides(1:end-1);

r = reshape(r, [prod(sides) basis_size]);
if ischar(gamma)
    if strcmp(gamma, 'even')
	gamma = 0;
    elseif strcmp(gamma, 'odd')
	gamma = inf;
    end
end

% It's more efficient in matlab code to do a small matrix multiplication
% "manually" in parallell over all the points than doing a multiple loop
% over the points and computing the matrix products "automatically". Thus we
% need to implement this differently for each tensor dimensionality.

switch basis_size
 case 3
  % 1D
  %
  % The tensor is of the form A*A'+gamma*b*b', where A and b are
  % composed from the elements in r according to
  % 
  % A=[3], b=[2].
  %
  % Thus (excluding gamma)
  %
  % T=[3*3+2*2].
  %
  % Comment: Orientation tensors in 1D are fairly pointless and only
  %          included here for completeness.

  T = zeros([prod(sides) 1 1]);
  
  if gamma == 0
      T(:,1,1) = r(:,3).^2;
  elseif gamma == inf
      T(:,1,1) = r(:,2).^2;
  else
      T(:,1,1) = r(:,3).^2 + gamma * r(:,2).^2;
  end
  T = reshape(T, [sides 1 1]);
  
 case 6
  % 2D
  %
  % The tensor is of the form A*A'+gamma*b*b', where A and b are
  % composed from the elements in r according to
  % 
  %   [4  6]    [2]
  % A=[6  5], b=[3]
  %
  % and the off-diagonal elements in A additionally are halved.
  %
  % Thus (excluding gamma)
  %
  %   [4*4+6*6+2*2 4*6+5*6+2*3]
  % T=[4*6+5*6+2*3 6*6+5*5+3*3]

  r(:,6) = r(:,6) / 2;
  T = zeros([prod(sides) 2 2]);
  
  if gamma == 0
      T(:,1,1) = r(:,4).^2 + r(:,6).^2;
      T(:,2,2) = r(:,5).^2 + r(:,6).^2;
      T(:,1,2) = (r(:,4) + r(:,5)) .* r(:,6);
  elseif gamma == inf
      T(:,1,1) = r(:,2).^2;
      T(:,2,2) = r(:,3).^2;
      T(:,1,2) = r(:,2) .* r(:,3);
  else
      T(:,1,1) = r(:,4).^2 + r(:,6).^2 + ...
		 gamma * r(:,2).^2;
      T(:,2,2) = r(:,5).^2 + r(:,6).^2 + ...
		 gamma * r(:,3).^2;
      T(:,1,2) = (r(:,4) + r(:,5)) .* r(:,6) + ...
		 gamma * r(:,2) .* r(:,3);
  end
  T(:,2,1) = T(:,1,2);
  T = reshape(T, [sides 2 2]);
  
 case 10
  % 3D
  %
  % The tensor is of the form A*A'+gamma*b*b', where A and b are
  % composed from the elements in r according to
  % 
  %   [5  8  9]    [2]
  % A=[8  6 10], b=[3]
  %   [9 10  7]    [4]
  %
  % and the off-diagonal elements in A additionally are halved.
  %
  % Thus (excluding gamma)
  %
  %   [5*5+8*8+9*9+2*2  5*8+6*8+9*10+2*3  5*9+8*10+7*9+2*4 ]
  % T=[5*8+6*8+9*10+2*3 8*8+6*6+10*10+3*3 8*9+6*10+7*10+3*4].
  %   [5*9+8*10+7*9+2*4 8*9+6*10+7*10+3*4 9*9+10*10+7*7+4*4]

  r(:,[8:10]) = r(:,[8:10]) / 2;
  T = zeros([prod(sides) 3 3]);
  
  if gamma == 0
      T(:,1,1) = r(:,5).^2 + r(:,8).^2 + r(:,9).^2;
      T(:,2,2) = r(:,8).^2 + r(:,6).^2 + r(:,10).^2;
      T(:,3,3) = r(:,9).^2 + r(:,10).^2 + r(:,7).^2;
      T(:,1,2) = (r(:,5) + r(:,6)) .* r(:,8) + r(:,9) .* r(:,10);
      T(:,1,3) = (r(:,5) + r(:,7)) .* r(:,9) + r(:,8) .* r(:,10);
      T(:,2,3) = (r(:,6) + r(:,7)) .* r(:,10) + r(:,8) .* r(:,9);
  elseif gamma == inf
      T(:,1,1) = r(:,2).^2;
      T(:,2,2) = r(:,3).^2;
      T(:,3,3) = r(:,4).^2;
      T(:,1,2) = r(:,2) .* r(:,3);
      T(:,1,3) = r(:,2) .* r(:,4);
      T(:,2,3) = r(:,3) .* r(:,4);
  else
      T(:,1,1) = r(:,5).^2 + r(:,8).^2 + r(:,9).^2 + ...
		 gamma * r(:,2).^2;
      T(:,2,2) = r(:,8).^2 + r(:,6).^2 + r(:,10).^2 + ...
		 gamma * r(:,3).^2;
      T(:,3,3) = r(:,9).^2 + r(:,10).^2 + r(:,7).^2 + ...
		 gamma * r(:,4).^2;
      T(:,1,2) = (r(:,5) + r(:,6)) .* r(:,8) + r(:,9) .* r(:,10) + ...
		 gamma * r(:,2) .* r(:,3);
      T(:,1,3) = (r(:,5) + r(:,7)) .* r(:,9) + r(:,8) .* r(:,10) + ...
		 gamma * r(:,2) .* r(:,4);
      T(:,2,3) = (r(:,6) + r(:,7)) .* r(:,10) + r(:,8) .* r(:,9) + ...
		 gamma * r(:,3) .* r(:,4);
  end
  T(:,2,1) = T(:,1,2);
  T(:,3,1) = T(:,1,3);
  T(:,3,2) = T(:,2,3);
  T = reshape(T, [sides 3 3]);
  
 case 15
  % 4D
  %
  % The tensor is of the form A*A'+gamma*b*b', where A and b are
  % composed from the elements in r according to
  % 
  %   [6  10 11 12]    [2]
  %   [10  7 13 14]    [3] 
  % A=[11 13  8 15], b=[4]
  %   [12 14 15  9]    [5]
  %
  % and the off-diagonal elements in A additionally are halved.
  %
  % Thus (excluding gamma)
  %
  %   [6*6+10*10+11*11+12*12+2*2 6*10+7*10+11*13+12*14+2*3
  %   [6*10+7*10+11*13+12*14+2*3 10*10+7*7+13*13+14*14+3*3
  % T=[6*11+10*13+8*11+12*15+2*4 10*11+7*13+8*13+14*15+3*4
  %   [6*12+10*14+11*15+9*12+2*5 10*12+7*14+13*15+9*14+3*5
  %
  %    6*11+10*13+8*11+12*15+2*4 6*12+10*14+11*15+9*12+2*5]
  %    10*11+7*13+8*13+14*15+3*4 10*12+7*14+13*15+9*14+3*5]
  %    11*11+13*13+8*8+15*15+4*4 11*12+13*14+8*15+9*15+4*5].
  %    11*12+13*14+8*15+9*15+4*5 12*12+14*14+15*15+9*9+5*5]

  r(:,[10:15]) = r(:,[10:15]) / 2;
  T = zeros([prod(sides) 4 4]);
  
  if gamma == 0
      T(:,1,1) = r(:,6).^2 + r(:,10).^2 + r(:,11).^2 + r(:,12).^2;
      T(:,2,2) = r(:,10).^2 + r(:,7).^2 + r(:,13).^2 + r(:,14).^2;
      T(:,3,3) = r(:,11).^2 + r(:,13).^2 + r(:,8).^2 + r(:,15).^2;
      T(:,4,4) = r(:,11).^2 + r(:,13).^2 + r(:,8).^2 + r(:,15).^2;
      T(:,1,2) = (r(:,6) + r(:,7)) .* r(:,10) + ...
		 r(:,11) .* r(:,13) + r(:,12) .* r(:,14);
      T(:,1,3) = (r(:,6) + r(:,8)) .* r(:,11) + ...
		 r(:,10) .* r(:,13) + r(:,12) .* r(:,15);
      T(:,1,4) = (r(:,6) + r(:,9)) .* r(:,12) + ...
		 r(:,10) .* r(:,14) + r(:,11) .* r(:,15);
      T(:,2,3) = (r(:,7) + r(:,8)) .* r(:,13) + ...
		 r(:,10) .* r(:,11) + r(:,14) .* r(:,15);
      T(:,2,4) = (r(:,7) + r(:,9)) .* r(:,14) + ...
		 r(:,10) .* r(:,12) + r(:,13) .* r(:,15);
      T(:,3,4) = (r(:,8) + r(:,9)) .* r(:,15) + ...
		 r(:,11) .* r(:,12) + r(:,13) .* r(:,14);
  elseif gamma == inf
      T(:,1,1) = r(:,2).^2;
      T(:,2,2) = r(:,3).^2;
      T(:,3,3) = r(:,4).^2;
      T(:,4,4) = r(:,4).^2;
      T(:,1,2) = r(:,2) .* r(:,3);
      T(:,1,3) = r(:,2) .* r(:,4);
      T(:,1,4) = r(:,2) .* r(:,5);
      T(:,2,3) = r(:,3) .* r(:,4);
      T(:,2,4) = r(:,3) .* r(:,5);
      T(:,3,4) = r(:,4) .* r(:,5);
  else
      T(:,1,1) = r(:,6).^2 + r(:,10).^2 + r(:,11).^2 + r(:,12).^2 + ...
		 gamma * r(:,2).^2;
      T(:,2,2) = r(:,10).^2 + r(:,7).^2 + r(:,13).^2 + r(:,14).^2 + ...
		 gamma * r(:,3).^2;
      T(:,3,3) = r(:,11).^2 + r(:,13).^2 + r(:,8).^2 + r(:,15).^2 + ...
		 gamma * r(:,4).^2;
      T(:,4,4) = r(:,11).^2 + r(:,13).^2 + r(:,8).^2 + r(:,15).^2 + ...
		 gamma * r(:,4).^2;
      T(:,1,2) = (r(:,6) + r(:,7)) .* r(:,10) + ...
		 r(:,11) .* r(:,13) + r(:,12) .* r(:,14) + ...
		 gamma * r(:,2) .* r(:,3);
      T(:,1,3) = (r(:,6) + r(:,8)) .* r(:,11) + ...
		 r(:,10) .* r(:,13) + r(:,12) .* r(:,15) + ...
		 gamma * r(:,2) .* r(:,4);
      T(:,1,4) = (r(:,6) + r(:,9)) .* r(:,12) + ...
		 r(:,10) .* r(:,14) + r(:,11) .* r(:,15) + ...
		 gamma * r(:,2) .* r(:,5);
      T(:,2,3) = (r(:,7) + r(:,8)) .* r(:,13) + ...
		 r(:,10) .* r(:,11) + r(:,14) .* r(:,15) + ...
		 gamma * r(:,3) .* r(:,4);
      T(:,2,4) = (r(:,7) + r(:,9)) .* r(:,14) + ...
		 r(:,10) .* r(:,12) + r(:,13) .* r(:,15) + ...
		 gamma * r(:,3) .* r(:,5);
      T(:,3,4) = (r(:,8) + r(:,9)) .* r(:,15) + ...
		 r(:,11) .* r(:,12) + r(:,13) .* r(:,14) + ...
		 gamma * r(:,4) .* r(:,5);
  end
  T(:,2,1) = T(:,1,2);
  T(:,3,1) = T(:,1,3);
  T(:,4,1) = T(:,1,4);
  T(:,3,2) = T(:,2,3);
  T(:,4,2) = T(:,2,4);
  T(:,4,3) = T(:,3,4);
  T = reshape(T, [sides 4 4]);

 otherwise
  error('Strange format for r.')
end

