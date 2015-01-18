function [sequence, mask, correct_flow] = read_yosemite
% [SEQUENCE, MASK, CORRECT_FLOW] = READ_YOSEMITE
%
% Read 15 frames from the Yosemite sequence. These are originally
% numbered yos2--yos16, with the correct flow being known for the middle
% frame yos9.
%
% SEQUENCE     - 252x316x15 array with gray values.
% MASK         - 252x316 logical array which is one everywhere except
%                for the sky.
% CORRECT_FLOW - The correct flow field for yos9, which is available
%                as SEQUENCE(:,:,8).
%
% The convention used for the motion is that the first coordinate is
% directed downwards in the image and the second coordinate is
% directed rightwards. This is consistent with the indexing into the
% image being stored as a matlab matrix.
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

s = which('read_yosemite');
slashes = findstr(s, '/');
yosemite_dir = [s(1:slashes(end)) 'yosemite_sequence/'];

% Read a sample frame and allocate space for 15 frames.
a = imread([yosemite_dir 'yos9.tif']);
sides = size(a);
sequence = zeros([sides 15]);

% Read the relevant part of the sequence.
for i=1:15
    sequence(:,:,i) = double(imread([yosemite_dir 'yos' int2str(i+1) '.tif']));
end

% Create the mask for the non-sky.
load([yosemite_dir 'correct_flow']);
mask = (correct_flow(:,:,1) ~= 2.0);

% Permute the directions to agree with our conventions.
correct_flow(:,:,[1 2]) = correct_flow(:,:,[2 1]);
correct_flow(:,:,1) = -correct_flow(:,:,1);
