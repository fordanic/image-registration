function test_yosemite
% TEST_YOSEMITE
%
% This function reproduces the results in the last five rows of table 6.2,
% page 102, of Gunnar Farnebäck's thesis "Polynomial Expansion for
% Orientation and Motion Estimation". The results in table 7.1, page 126,
% are also reproduced.
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

% Read the Yosemite sequence.
[sequence, mask, correct_flow] = read_yosemite;

disp('Yosemite, fast algorithm, constant motion:')
test_velocity1(sequence, mask, correct_flow, 9, 1.4, 1/32, ...
	       'constant', 15, 3.5);

disp('Yosemite, fast algorithm, affine motion:')
test_velocity1(sequence, mask, correct_flow, 11, 1.6, 1/256, ...
	       'affine', 41, 6.5);

% This is not included in table 6.2
disp('Yosemite, segmentation algorithm, affine motion, single run:')
test_velocity2(sequence, mask, correct_flow, 9, 1.4, 1/8, ...
	       'affine', 500, 0.06);

disp('Yosemite, segmentation algorithm, affine motion, 11 runs averaged:')
test_velocity2(sequence, mask, correct_flow, 9, 1.4, 1/8, ...
	       'affine', [400:20:600], 0.06);

disp('Yosemite, displacement algorithm, constant motion, 1 iteration:')
test_velocity3(sequence, mask, correct_flow, 11, 41, 'constant', 'accurate');

disp('Yosemite, displacement algorithm, constant motion, 3 iterations:')
test_velocity3(sequence, mask, correct_flow, [11 11 11], [41 41 41], ...
	       'constant', 'accurate');

% This result is slightly better than the one in table 7.1. The difference
% is that a slightly larger averaging kernel is used here, although with the
% same sigma.
disp('Yosemite, displacement algorithm, affine motion, 1 iteration:')
test_velocity3(sequence, mask, correct_flow, 11, 41, 'affine', 'accurate');

disp('Yosemite, displacement algorithm, affine motion, 3 iterations:')
test_velocity3(sequence, mask, correct_flow, [11 11 11], [41 41 41], ...
	       'affine', 'accurate');

% These are the same as the previous four except that the ordering of the
% two frames is different. These results are equal or better than the ones
% in table 7.1.
disp('Yosemite, displacement algorithm, constant motion, 1 iteration:')
test_velocity4(sequence, mask, correct_flow, 11, 41, 'constant', 'accurate');

disp('Yosemite, displacement algorithm, constant motion, 3 iterations:')
test_velocity4(sequence, mask, correct_flow, [11 11 11], [41 41 41], ...
	       'constant', 'accurate');

disp('Yosemite, displacement algorithm, affine motion, 1 iteration:')
test_velocity4(sequence, mask, correct_flow, 11, 41, 'affine', 'accurate');

disp('Yosemite, displacement algorithm, affine motion, 3 iterations:')
test_velocity4(sequence, mask, correct_flow, [11 11 11], [41 41 41], ...
	       'affine', 'accurate');
return


% Fast velocity algorithm.
function test_velocity1(sequence, mask, correct_flow, kernelsize, sigma, ...
			gamma, model, kernelsize_avg, sigma_avg)

% Extract the size of the frames.
sides = size(sequence);
middle = (sides(3)+1)/2;
sides = sides(1:2);

% Region of interest for tensor computation. We only compute the tensor
% field for the middle frame.
roi = [[[1;1] sides'];[middle middle]];

% Compute tensors.
options = struct('sigma', sigma, 'gamma', gamma);
T = make_tensors_fast(sequence, kernelsize, roi, options);

% Certainty mask for averaging. Close to the border the tensors will be
% affected by edge effects, so we ignore them.
bwidth = (kernelsize-1)/2;
mask2 = ones(sides);
mask2([1:bwidth,end-bwidth+1:end],:) = 0;
mask2(:,[1:bwidth,end-bwidth+1:end]) = 0;

% Compute velocity from a tensor field using the fast algorithm and the
% specified motion model.
[v, c] = velocity_from_tensors(T, model, kernelsize_avg, sigma_avg, mask2);

% Evaluate the estimated velocity.
evaluate_velocity(v, correct_flow, mask, c);
return


% Segmentation velocity algorithm.
function test_velocity2(sequence, mask, correct_flow, kernelsize, sigma, ...
			gamma, model, cand_sizes, lambda)

% Extract the size of the frames.
sides = size(sequence);
middle = (sides(3)+1)/2;
sides = sides(1:2);

% Region of interest for tensor computation. We only compute the tensor
% field for the middle frame.
roi1=[[[1;1] sides'];[middle middle]];

% Close to the edges we improve the tensor estimates by using the separable
% normalized convolution method for polynomial expansion. We only do this
% close to the edges because it is slower and gives the same results for the
% central part of the signal.
bwidth = (kernelsize-1)/2;
roi2 = [1 bwidth;1 sides(2);middle middle];
roi3 = [sides(1)-bwidth+1 sides(1);1 sides(2);middle middle];
roi4 = [bwidth+1 sides(1)-bwidth;1 bwidth;middle middle];
roi5 = [bwidth+1 sides(1)-bwidth;sides(2)-bwidth+1 sides(2);middle middle];

options = struct('sigma',sigma,'gamma',gamma);

% Compute tensors.
T = make_tensors_fast(sequence, kernelsize, roi1, options);

cert = ones(size(sequence));

clear options;
options.region_of_interest = roi2;
r = polyexp(sequence, cert, 'quadratic', kernelsize, sigma, options);
T(1:bwidth,:,:,:) = poly_to_tensor(r, gamma);

options.region_of_interest = roi3;
r = polyexp(sequence, cert, 'quadratic', kernelsize, sigma, options);
T(end-bwidth+1:end,:,:,:) = poly_to_tensor(r, gamma);

options.region_of_interest = roi4;
r = polyexp(sequence, cert, 'quadratic', kernelsize, sigma, options);
T(bwidth+1:end-bwidth,1:bwidth,:,:) = poly_to_tensor(r, gamma);

options.region_of_interest = roi5;
r = polyexp(sequence, cert, 'quadratic', kernelsize, sigma, options);
T(bwidth+1:end-bwidth,end-bwidth+1:end,:,:) = poly_to_tensor(r, gamma);

% Remove the isotropic part of the tensor.
T = remove_isotropic(T);

if strcmp(model, 'constant')
    model = 1;
elseif strcmp(model, 'affine')
    model = 2;
else
    model = 4; % Eightparam
end

v = zeros([sides 2]);
for k = 1:length(cand_sizes)
    options = struct('lambda', lambda, 'motion_models', model);
    options.minsize = cand_sizes(k);
    options.maxsize = cand_sizes(k);
    this_v = velocity_segment(T, options);
    v = v + this_v;
end
v = v / length(cand_sizes);

% Evaluate the estimated velocity.
evaluate_velocity(v, correct_flow, mask, ones(size(v)), 1);
return


% Displacement estimation algorithm.
function test_velocity3(sequence, mask, correct_flow, kernelsizes1, ...
			kernelsizes2, model, method)

% Extract the size of the frames.
sides = size(sequence);
middle = (sides(3)+1)/2;
sides = sides(1:2);

% Compute displacement from the middle frame and the one before.
% Notice that section 7.9.2 claims that the middle frame and the *following*
% frame are used. That is incorrect and was caused by a bug in the code
% (an embarrasing off by one error).
im1 = sequence(:,:,middle-1);
im2 = sequence(:,:,middle);

% Estimate the displacement.
[v, c] = estimate_displacement(im1, im2, kernelsizes1, kernelsizes2, ...
			       model, method);

% Evaluate the estimated velocity.
evaluate_velocity(v, correct_flow, mask, c, [1 0.9 0.7 0.5 0.3]);
return


% Displacement estimation algorithm.
function test_velocity4(sequence, mask, correct_flow, kernelsizes1, ...
			kernelsizes2, model, method)

% Extract the size of the frames.
sides = size(sequence);
middle = (sides(3)+1)/2;
sides = sides(1:2);

% Compute displacement from the middle frame and the one before.
im1 = sequence(:,:,middle-1);
im2 = sequence(:,:,middle);

% Estimate the displacement. Compared to test_velocity3 we reverse the order
% of the frames. This makes a subtle difference when we use iterated
% estimation.
[v, c] = estimate_displacement(im2, im1, kernelsizes1, kernelsizes2, ...
			       model, method);

% Evaluate the estimated velocity. Since the frame order was reversed we
% need to negate the estimated displacement field.
evaluate_velocity(-v, correct_flow, mask, c, [1 0.9 0.7 0.5 0.3]);
return


