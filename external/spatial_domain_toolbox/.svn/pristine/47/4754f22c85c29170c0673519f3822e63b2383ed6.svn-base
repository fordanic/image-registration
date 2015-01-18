function [v, cout] = velocity_from_tensors(T, model, kernelsize, sigma, cin)
% [V, COUT] = VELOCITY_FROM_TENSORS(T, MODEL, KERNELSIZE, SIGMA, CIN)
%
% Compute a velocity field from a tensor field, using the algorithm
% described in section 6.5 of Gunnar Farnebäck's thesis "Polynomial
% Expansion for Orientation and Motion Estimation". Normalized averaging is
% performed separably with a Gaussian filter.
%
% T          - Tensor field of size height x width x 3 x 3.
% MODEL      - Velocity field model. Supported are 'constant', 'affine',
%              and 'eightparam', with definitions according to equations
%              6.4, 6.5, and 6.6.
% KERNELSIZE - Spatial size of the Gaussian averaging filter.
% SIGMA      - Standard deviation of the Gaussian averaging filter
% CIN        - Input certainty for the tensor field, used in the normalized
%              averaging step. Typically cin should be reduced close to the
%              borders of the tensor field.
% V          - Estimated velocity field of size height x width x 2
% COUT       - Output confidence field, according to equation 6.23. Notice
%              that a small value indicates a high confidence.
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

% Extract the spatial size of the tensor field
sides = size(T);
sides = sides(1:2);

% Remove the isotropic part of the tensors.
T = remove_isotropic(T);

switch model
  case 'constant'
    % Tensor averaging, which in this case is equivalent with computation of Q.
    T = T(:,:,[1 5 9 2 3 6]);
    app = gaussian_app(kernelsize, 1, sigma);
    T = conv3(conv3(T.*repmat(cin, [1 1 6]), app), app') ./ ...
	(eps + repmat(conv3(conv3(cin, app), app'), [1 1 6])); 
    T = reshape(T(:,:,[1 4 5;4 2 6;5 6 3]), [sides 3 3]);

    %Solve the equation Qv=-q.
    a = T(:,:,1,1);
    b = T(:,:,2,2);
    c = T(:,:,1,2);
    d = T(:,:,1,3);
    e = T(:,:,2,3);
    % Q = [a c;c b], inv(Q) = [b -c;-c a]/(a*b-c^2), q = [d;e]
    v = -cat(3, d.*b-c.*e, a.*e-c.*d) ./ (eps + repmat(a.*b-c.^2, [1 1 2]));

    %Compute output certainty
    if nargout > 1
	cout = T(:,:,3,3) + d.*v(:,:,1) + e.*v(:,:,2);
    end
 
  case 'affine'
    % Build the quadratic forms for the affine motion
    [x,y] = meshgrid(1:sides(2),1:sides(1));
    Q = zeros([sides 25]);
    Q(:,:,1) = x.*x.*T(:,:,1,1);
    Q(:,:,2) = y.*y.*T(:,:,1,1);
    Q(:,:,3) = T(:,:,1,1);
    Q(:,:,4) = x.*x.*T(:,:,2,2);
    Q(:,:,5) = y.*y.*T(:,:,2,2);
    Q(:,:,6) = T(:,:,2,2);
    Q(:,:,7) = T(:,:,3,3);
    
    Q(:,:,8) = x.*y.*T(:,:,1,1);
    Q(:,:,9) = x.*T(:,:,1,1);
    Q(:,:,10) = x.*x.*T(:,:,1,2);
    Q(:,:,11) = x.*y.*T(:,:,1,2);
    Q(:,:,12) = x.*T(:,:,1,2);
    Q(:,:,13) = x.*T(:,:,1,3);
    
    Q(:,:,14) = y.*T(:,:,1,1);
    Q(:,:,15) = y.*y.*T(:,:,1,2);
    Q(:,:,16) = y.*T(:,:,1,2);
    Q(:,:,17) = y.*T(:,:,1,3);
    
    Q(:,:,18) = T(:,:,1,2);
    Q(:,:,19) = T(:,:,1,3);
    
    Q(:,:,20) = x.*y.*T(:,:,2,2);
    Q(:,:,21) = x.*T(:,:,2,2);
    Q(:,:,22) = x.*T(:,:,2,3);
    
    Q(:,:,23) = y.*T(:,:,2,2);
    Q(:,:,24) = y.*T(:,:,2,3);
    
    Q(:,:,25) = T(:,:,2,3);

    app = gaussian_app(kernelsize, 1, sigma);
    Q = conv3(conv3(Q.*repmat(cin, [1 1 25]), app), app') ./ ...
	(eps + repmat(conv3(conv3(cin, app), app'), [1 1 25])); 
    Q = reshape(Q(:,:,[ 1  8  9 10 11 12 13;
			8  2 14 11 15 16 17;
			9 14  3 12 16 18 19;
			10 11 12  4 20 21 22;
			11 15 16 20  5 23 24;
			12 16 18 21 23  6 25;
			13 17 19 22 24 25  7]),[sides 7 7]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The code above is somewhat obfuscated in order to improve speed. What
    % we really want to do is to compute Q = S^T T S as in equation 5.17 and
    % Q_{tot} according to equation 5.31, with the affine motion model.
    % A more direct, but slower, implementation is shown below:
    %
    % [x,y] = meshgrid(1:sides(2),1:sides(1));
    % Q = zeros([sides 7 7]);
    % Q(:,:,1,1) = x.*x.*T(:,:,1,1);
    % Q(:,:,2,2) = y.*y.*T(:,:,1,1);
    % Q(:,:,3,3) = T(:,:,1,1);
    % Q(:,:,4,4) = x.*x.*T(:,:,2,2);
    % Q(:,:,5,5) = y.*y.*T(:,:,2,2);
    % Q(:,:,6,6) = T(:,:,2,2);
    % Q(:,:,7,7) = T(:,:,3,3);
    % 
    % Q(:,:,1,2) = x.*y.*T(:,:,1,1);
    % Q(:,:,1,3) = x.*T(:,:,1,1);
    % Q(:,:,1,4) = x.*x.*T(:,:,1,2);
    % Q(:,:,1,5) = x.*y.*T(:,:,1,2);
    % Q(:,:,1,6) = x.*T(:,:,1,2);
    % Q(:,:,1,7) = x.*T(:,:,1,3);
    % 
    % Q(:,:,2,3) = y.*T(:,:,1,1);
    % Q(:,:,2,4) = x.*y.*T(:,:,1,2);
    % Q(:,:,2,5) = y.*y.*T(:,:,1,2);
    % Q(:,:,2,6) = y.*T(:,:,1,2);
    % Q(:,:,2,7) = y.*T(:,:,1,3);
    % 
    % Q(:,:,3,4) = x.*T(:,:,1,2);
    % Q(:,:,3,5) = y.*T(:,:,1,2);
    % Q(:,:,3,6) = T(:,:,1,2);
    % Q(:,:,3,7) = T(:,:,1,3);
    % 
    % Q(:,:,4,5) = x.*y.*T(:,:,2,2);
    % Q(:,:,4,6) = x.*T(:,:,2,2);
    % Q(:,:,4,7) = x.*T(:,:,2,3);
    % 
    % Q(:,:,5,6) = y.*T(:,:,2,2);
    % Q(:,:,5,7) = y.*T(:,:,2,3);
    % 
    % Q(:,:,6,7) = T(:,:,2,3);
    % 
    % Q(:,:,2,1) = Q(:,:,1,2);
    % Q(:,:,3,1) = Q(:,:,1,3);
    % Q(:,:,4,1) = Q(:,:,1,4);
    % Q(:,:,5,1) = Q(:,:,1,5);
    % Q(:,:,6,1) = Q(:,:,1,6);
    % Q(:,:,7,1) = Q(:,:,1,7);
    % Q(:,:,3,2) = Q(:,:,2,3);
    % Q(:,:,4,2) = Q(:,:,2,4);
    % Q(:,:,5,2) = Q(:,:,2,5);
    % Q(:,:,6,2) = Q(:,:,2,6);
    % Q(:,:,7,2) = Q(:,:,2,7);
    % Q(:,:,4,3) = Q(:,:,3,4);
    % Q(:,:,5,3) = Q(:,:,3,5);
    % Q(:,:,6,3) = Q(:,:,3,6);
    % Q(:,:,7,3) = Q(:,:,3,7);
    % Q(:,:,5,4) = Q(:,:,4,5);
    % Q(:,:,6,4) = Q(:,:,4,6);
    % Q(:,:,7,4) = Q(:,:,4,7);
    % Q(:,:,6,5) = Q(:,:,5,6);
    % Q(:,:,7,5) = Q(:,:,5,7);
    % Q(:,:,7,6) = Q(:,:,6,7);
    % 
    % Q = Q(:,:,:);
    % app = gaussian_app(kernelsize,1,sigma);
    % Q = conv3(conv3(Q.*repmat(cin,[1 1 49]),app),app') ./ ...
    %     (eps+repmat(conv3(conv3(cin,app),app'),[1 1 49])); 
    % Q = reshape(Q,[sides 7 7]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Solve the equation Qp=-q.
    p = min_quadform(Q);
    
    v = zeros([sides 2]);
    v(:,:,1) = sum(p(:,:,1:3).*cat(3,x,y,ones(sides)),3);
    v(:,:,2) = sum(p(:,:,4:6).*cat(3,x,y,ones(sides)),3);
    if nargout > 1
    	cout = Q(:,:,7,7) + sum(p.*Q(:,:,[1:6],7),3);
    end
    
  case 'eightparam'
    % Build the quadratic forms for the eight parameter motion
    [x,y] = meshgrid(1:sides(2), 1:sides(1));
    Q = zeros([sides 9 9]);

    Q(:,:,1,1) = T(:,:,1,1);
    Q(:,:,2,2) = x.^2.*T(:,:,1,1);
    Q(:,:,3,3) = y.^2.*T(:,:,1,1);
    Q(:,:,4,4) = T(:,:,2,2);
    Q(:,:,5,5) = x.^2.*T(:,:,2,2);
    Q(:,:,6,6) = y.^2.*T(:,:,2,2);
    Q(:,:,7,7) = x.^4.*T(:,:,1,1) + 2*x.^3.*y.*T(:,:,1,2) + ...
		 x.^2.*y.^2.*T(:,:,2,2);
    Q(:,:,8,8) = x.^2.*y.^2.*T(:,:,1,1) + 2*x.*y.^3.*T(:,:,1,2) + ...
		 y.^4.*T(:,:,2,2);
    Q(:,:,9,9) = T(:,:,3,3);
    	  
    Q(:,:,1,2) = x.*T(:,:,1,1);
    Q(:,:,1,3) = y.*T(:,:,1,1);
    Q(:,:,1,4) = T(:,:,1,2);
    Q(:,:,1,5) = x.*T(:,:,1,2);
    Q(:,:,1,6) = y.*T(:,:,1,2);
    Q(:,:,1,7) = x.^2.*T(:,:,1,1) + x.*y.*T(:,:,1,2);
    Q(:,:,1,8) = x.*y.*T(:,:,1,1) + y.^2.*T(:,:,1,2);
    Q(:,:,1,9) = T(:,:,1,3);
    Q(:,:,2,3) = x.*y.*T(:,:,1,1);
    Q(:,:,2,4) = x.*T(:,:,1,2);
    Q(:,:,2,5) = x.^2.*T(:,:,1,2);
    Q(:,:,2,6) = x.*y.*T(:,:,1,2);
    Q(:,:,2,7) = x.^3.*T(:,:,1,1) + x.^2.*y.*T(:,:,1,2);
    Q(:,:,2,8) = x.^2.*y.*T(:,:,1,1) + x.*y.^2.*T(:,:,1,2);
    Q(:,:,2,9) = x.*T(:,:,1,3);
    Q(:,:,3,4) = y.*T(:,:,1,2);
    Q(:,:,3,5) = x.*y.*T(:,:,1,2);
    Q(:,:,3,6) = y.^2.*T(:,:,1,2);
    Q(:,:,3,7) = x.^2.*y.*T(:,:,1,1) + x.*y.^2.*T(:,:,1,2);
    Q(:,:,3,8) = x.*y.^2.*T(:,:,1,1) + y.^3.*T(:,:,1,2);
    Q(:,:,3,9) = y.*T(:,:,1,3);
    Q(:,:,4,5) = x.*T(:,:,2,2);
    Q(:,:,4,6) = y.*T(:,:,2,2);
    Q(:,:,4,7) = x.^2.*T(:,:,1,2) + x.*y.*T(:,:,2,2);
    Q(:,:,4,8) = x.*y.*T(:,:,1,2) + y.^2.*T(:,:,2,2);
    Q(:,:,4,9) = T(:,:,2,3);
    Q(:,:,5,6) = x.*y.*T(:,:,2,2);
    Q(:,:,5,7) = x.^3.*T(:,:,1,2) + x.^2.*y.*T(:,:,2,2);
    Q(:,:,5,8) = x.^2.*y.*T(:,:,1,2) + x.*y.^2.*T(:,:,2,2);
    Q(:,:,5,9) = x.*T(:,:,2,3);
    Q(:,:,6,7) = x.^2.*y.*T(:,:,1,2) + x.*y.^2.*T(:,:,2,2);
    Q(:,:,6,8) = x.*y.^2.*T(:,:,1,2) + y.^3.*T(:,:,2,2);
    Q(:,:,6,9) = y.*T(:,:,2,3);
    Q(:,:,7,8) = x.^3.*y.*T(:,:,1,1) + 2*x.^2.*y.^2.*T(:,:,1,2) + ...
		 x.*y.^3.*T(:,:,2,2);
    Q(:,:,7,9) = x.^2.*T(:,:,1,3) + x.*y.*T(:,:,2,3);
    Q(:,:,8,9) = x.*y.*T(:,:,1,3) + y.^2.*T(:,:,2,3);
    
    Q(:,:,2,1) = Q(:,:,1,2);
    Q(:,:,3,1) = Q(:,:,1,3);
    Q(:,:,4,1) = Q(:,:,1,4);
    Q(:,:,5,1) = Q(:,:,1,5);
    Q(:,:,6,1) = Q(:,:,1,6);
    Q(:,:,7,1) = Q(:,:,1,7);
    Q(:,:,8,1) = Q(:,:,1,8);
    Q(:,:,9,1) = Q(:,:,1,9);
    Q(:,:,3,2) = Q(:,:,2,3);
    Q(:,:,4,2) = Q(:,:,2,4);
    Q(:,:,5,2) = Q(:,:,2,5);
    Q(:,:,6,2) = Q(:,:,2,6);
    Q(:,:,7,2) = Q(:,:,2,7);
    Q(:,:,8,2) = Q(:,:,2,8);
    Q(:,:,9,2) = Q(:,:,2,9);
    Q(:,:,4,3) = Q(:,:,3,4);
    Q(:,:,5,3) = Q(:,:,3,5);
    Q(:,:,6,3) = Q(:,:,3,6);
    Q(:,:,7,3) = Q(:,:,3,7);
    Q(:,:,8,3) = Q(:,:,3,8);
    Q(:,:,9,3) = Q(:,:,3,9);
    Q(:,:,5,4) = Q(:,:,4,5);
    Q(:,:,6,4) = Q(:,:,4,6);
    Q(:,:,7,4) = Q(:,:,4,7);
    Q(:,:,8,4) = Q(:,:,4,8);
    Q(:,:,9,4) = Q(:,:,4,9);
    Q(:,:,6,5) = Q(:,:,5,6);
    Q(:,:,7,5) = Q(:,:,5,7);
    Q(:,:,8,5) = Q(:,:,5,8);
    Q(:,:,9,5) = Q(:,:,5,9);
    Q(:,:,7,6) = Q(:,:,6,7);
    Q(:,:,8,6) = Q(:,:,6,8);
    Q(:,:,9,6) = Q(:,:,6,9);
    Q(:,:,8,7) = Q(:,:,7,8);
    Q(:,:,9,7) = Q(:,:,7,9);
    Q(:,:,9,8) = Q(:,:,8,9);
    
    Q = Q(:,:,:);
    app = gaussian_app(kernelsize, 1, sigma);
    Q = conv3(conv3(Q.*repmat(cin, [1 1 81]), app), app') ./ ...
	(eps + repmat(conv3(conv3(cin, app), app'), [1 1 81]));
    Q = reshape(Q, [sides 9 9]);
    
    % Solve the equation Qp = -q.
    p = min_quadform(Q);
    v = zeros([sides 2]);
    v(:,:,1) = sum(p(:,:,[1 2 3 7 8]).*cat(3,ones(sides),x,y,x.^2,x.*y),3);
    v(:,:,2) = sum(p(:,:,4:8).*cat(3,ones(sides),x,y,x.*y,y.^2),3);
    if nargout > 1
	cout = Q(:,:,9,9) + sum(p.*Q(:,:,[1:8],9), 3);
    end
    
  otherwise
    error('invalid motion model')
end
