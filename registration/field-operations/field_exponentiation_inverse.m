function field = field_exponentiation_inverse(field)
% function field = field_exponentiation_inverse(field)
%
% Estimates the inverse field using a field exponentiation of
% the provided field.
%
% INPUT ARGUMENTS
% field 		- Field to perform inverse field exponentiation on
%
% OPTIONAL INPUT ARGUMENTS
% N/A
%
% OUTPUT ARGUMENTS
% field 		- Inverse field exponentiation

% Copyright (c) 2012 Daniel Forsberg
% danne.forsberg@outlook.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

fieldSqr = field{1}.^2;
for  k = 2 : length(field)
    fieldSqr = fieldSqr+field{k}.^2;
end

N = ceil(2 + log2(max(max(max(sqrt(fieldSqr)))))/2) +1;
if(N<1)
    N=1;
end

for k = 1 : length(field)
    field{k} = -field{k}*2^(-N);
end
for l = 1 : N
    for k = 1 : length(field)
        newField{k} = deformation(field{k},field,'linear');
        newField{k} = newField{k}+field{k};
    end
    field = newField;
end