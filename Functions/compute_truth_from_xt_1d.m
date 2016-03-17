function [ytruth_on_mcmcgrid] = compute_truth_from_xt_1d(x_grid, shape)
% Evaluate true convex function on x_grid (from compute_mesh_info_1d).
% NOTE: Functions appear in BOTH this file and make_noisy_convex_1d.m.
%
% Args:
%   x_grid: Array of points.
%   shape: String indicating which convex function to use.
%
% Returns:
%   ytruth_on_mcmcgrid: Matrix (nx1) of true convex response values.

n = length(x_grid);
ytruth_on_mcmcgrid = zeros(n, 1);

if strcmp(shape, 'parabola')
    for ii = 1:n
        ytruth_on_mcmcgrid(ii) = 1e-3 * x_grid(ii)^4;            
    end
    
else
    error('Shape should be either "parabola" or "parabola".')
end

end

