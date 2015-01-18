function z = outerprod(x, y, varargin)
% OUTERPROD - Multiply two or more arrays to form an outer product.
%             Singular dimensions are ignored and automatically squeezed.
%
% Author: Gunnar Farnebäck
%         Medical Informatics
%         Linköping University, Sweden
%         gf@isy.liu.se

dims = [size(x) size(y)];
dims = dims(dims ~= 1);
z = reshape(x(:) * y(:)', dims);

if length(varargin) > 0
    z = outerprod(z, varargin{:});
end
