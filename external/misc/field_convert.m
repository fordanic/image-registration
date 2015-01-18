function field = field_convert(inputField,dims)
% FIELD_CONVERT Converts cell field to matrix field and vice versa
%
% field = field_convert(inputField,dims)
%
% inputField    - Field to convert
% dims          - Number of dimensions
% field         - Converted field

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

if(iscell(inputField))
    if(nargin<2)
        dims = length(inputField);
    end

    if(dims==2)

        field = zeros([2,size(cell2mat(inputField(1)))],'double');
        field(1,:,:) = cell2mat(inputField(2));
        field(2,:,:) = cell2mat(inputField(1));

    elseif(dims==3)

        field = zeros([3,size(cell2mat(inputField(1)))],'double');
        field(1,:,:,:) = cell2mat(inputField(2));
        field(2,:,:,:) = cell2mat(inputField(1));
        field(3,:,:,:) = cell2mat(inputField(3));

    end

else
    field = cell(0);
    if(nargin<2)
        dims = ndims(inputField)-1;
    end

    if(dims==2)

        field{1} = squeeze(inputField(2,:,:));
        field{2} = squeeze(inputField(1,:,:));

    elseif(dims==3)

        field{1} = squeeze(inputField(2,:,:,:));
        field{2} = squeeze(inputField(1,:,:,:));
        field{3} = squeeze(inputField(3,:,:,:));

    end


end