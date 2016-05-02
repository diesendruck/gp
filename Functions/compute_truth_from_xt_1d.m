function [ytruth_on_grid] = compute_truth_from_xt_1d(x_grid, shape)
% Evaluate true convex function on x_grid (from compute_mesh_info_1d).
% NOTE: Functions appear in BOTH this file and make_noisy_convex_1d.m.
%
% Args:
%   x_grid: Array of points.
%   shape: String indicating which convex function to use.
%
% Returns:
%   ytruth_on_grid: Matrix (nx1) of true convex response values.

n = length(x_grid);
ytruth_on_grid = zeros(n, 1);

if strcmp(shape, '1d_parabola')
    for ii = 1:n
        ytruth_on_grid(ii) = 0.004 * (x_grid(ii)-7)^4;  
    end

elseif strcmp(shape, '1d_exponential')
    for ii = 1:n
        %ytruth_on_mcmcgrid(ii) = exp(x_grid(ii)) + 0.5*x_grid(ii)^2;            
        ytruth_on_grid(ii) = 4.5e-4 * exp(x_grid(ii));
    end
    
elseif strcmp(shape, '1d_negative_entropy')
    for ii = 1:n
        ytruth_on_grid(ii) = 0.4 * x_grid(ii) * log(x_grid(ii)+1); 
    end
    
else
    error('Shape not recognized.')
end

end

