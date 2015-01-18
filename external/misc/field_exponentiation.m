function field = field_exponentiation(field)
% FIELD_EXPONENTIATION Performs a field exponentiation on the provided field.
%
% field = field_exponentiation(field)
%
% field         - Field to exponentiate

%{ 
Copyright © 2011 by Université catholique de Louvain.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files, to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

REGGUI SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
%}

fieldSqr = field{1}.^2;
for  k = 2 : length(field)
    fieldSqr = fieldSqr+field{k}.^2;
end

N = ceil(2 + log2(max(sqrt(fieldSqr(:))))/2) + 1;
if(N<1)
    N=1;
end

for k = 1 : length(field)
    field{k} = field{k}*2^(-N);
end
for l = 1 : N
    for k = 1 : length(field)
        newField{k} = deformation(field{k},field,'linear');
        newField{k} = newField{k}+field{k};
    end
    field = newField;
end