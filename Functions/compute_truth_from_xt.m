function [ytruth_on_mcmcgrid] = compute_truth_from_xt(xt, shape)
% Evaluate true convex function on xt grid (from run_gpmc).

n = length(xt);
ytruth_on_mcmcgrid = zeros(n, 1);

if strcmp(shape, 'trough')  % Convex wrt first dim of x. Other variable is random depth.
    for ii = 1:n
        ytruth_on_mcmcgrid(ii) = 0.25 * xt(ii, 1)^2;            
    end
    
elseif strcmp(shape, 'paraboloid')
    for jj = 1:n
        ytruth_on_mcmcgrid(jj) = norm(xt(jj, :))^2;
    end
else
    error('Shape should be either "trough" or "paraboloid".')
end

end

