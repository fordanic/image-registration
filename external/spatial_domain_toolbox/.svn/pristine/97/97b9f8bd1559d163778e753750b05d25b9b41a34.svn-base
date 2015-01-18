function cout = polyexp_cout_helper(r, G, options)
% Helper function for polyexp.m.
% This would have been a local function in polyexp.m if that had worked
% in Matlab (tested with 5.3 and 6.5).
%
% Author: Gunnar Farnebäck
%         Computer Vision Laboratory
%         Linköping University, Sweden
%         gf@isy.liu.se

h = G * r;

if (isfield(options, 'cout_data'))
	cout = feval(options.cout_func, G, G, h, r, options.cout_data);
else
	cout = feval(options.cout_func, G, G, h, r);
end
