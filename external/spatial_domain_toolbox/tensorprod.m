function [t, outind] = tensorprod(varargin)
% TENSORPROD - Multiply two or more arrays while summing over repeated
%              indices. The indices for each array are indicated by strings.
%
% For example,
%
% C = tensorprod(A, 'ij', B, 'jk');
%
% computes a normal matrix product. The names of the indices are limited to
% a single character but are of no importance except to indicate repeated
% indices. The output indices come in the order they appear in the input
% arrays. Thus
%
% C = tensorprod(A, 'ki', B, 'ij');
%
% is exactly the same as the previous example. The output indices can be
% obtained by
%
% [C, ind] = tensorprod(A, 'ki', B, 'ij');
%
% and would be 'kj' in the latter example, whereas it would be 'ik' in the
% first example.
%
% No memory is spent on possible intermediate results but the speed is not
% as fast as a special-written implementation could be.
%
% BUGS: Imaginary parts of arrays are silently ignored.
% 
% Author: Gunnar Farnebäck
%         Medical Informatics
%         Linköping University, Sweden
%         gf@isy.liu.se

% Check for suspected call of old definition.
string_found = 0;
for k = 1:length(varargin)
    if ischar(varargin{k})
        string_found = 1;
    end
end
if ~string_found
    error(sprintf(['tensorprod has been changed. Update your code for ' ...
                   'the new definition\nor convert to the new function' ...
                   ' outerprod.']))
end

if mod(length(varargin), 2) == 1
    error('There must be an even number of arguments.');
end

translation_table = zeros(256, 1);
index_number = 1;
output_indices = [];
index_sizes = [];
arrays = cell(length(varargin) / 2, 1);
index_vectors = cell(length(varargin) / 2, 1);
index_vectors = {};

for k = 1:2:length(varargin)
    x = varargin{k};
    ind = varargin{k + 1};
    if ~ischar(ind)
        error(sprintf('Argument %d expected to be a string.', k + 1))
    end
    
    n = ndims(x);
    if n == 2 && size(x, 2) == 1
        n = 1;
    end
    if n ~= length(ind)
        error(sprintf(['Mismatch between array dimension and number of ' ...
                       'indices, args %d and %d.'], k, k + 1))
    end
    
    indices = zeros(n, 1);
    for m = 1:n
        c = double(ind(m));
        v = translation_table(c);
        if v == 0
            v = index_number;
            index_number = index_number + 1;
            translation_table(c) = v;
            output_indices = [output_indices v];
            index_sizes(v) = size(x, m);
        else
            output_indices = output_indices(output_indices ~= v);
            if (index_sizes(v) ~= size(x, m))
                error('All repeated indices must have the same length.');
            end
        end
        indices(m) = v;
    end
    
    arrays{(k + 1) / 2} = x;
    index_vectors{(k + 1) / 2} = indices;
end

if nargout > 1
    outind = '';
    for k = 1:length(output_indices)
        outind = [outind char(find(translation_table == output_indices(k)))];
    end
end

% If you need to speed up repeated tensorprod computations you can call
% tensorprodc directly with the arguments assembled here, but beware that it
% doesn't do any kind of error checking and might crash the whole Matlab
% process if you mess it up.
t = tensorprodc(index_sizes, output_indices, arrays{:}, index_vectors{:});
