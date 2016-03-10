function [ytruth_on_mcmcgrid] = compute_truth_from_xt(xt, shape)
% Evaluate true convex function on xt grid (from run_gpmc).
%
% Args:
%   xt: Meshgrid matrix of points.
%   shape: String indicating which convex function to use.
%
% Returns:
%   ytruth_on_mcmcgrid: Matrix (nx1) of true convex response values.

n = length(xt);
ytruth_on_mcmcgrid = zeros(n, 1);

if strcmp(shape, 'trough')  % Convex wrt first dim of x. Other variable is random depth.
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = 0.25*x1^2;            
    end
    
elseif strcmp(shape, 'paraboloid')
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = 0.05*sum(x1^2 + x2^2);
    end
    
elseif strcmp(shape, 'hand')
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = 0.1*x1^2 - log(x2) + 10;
    end
    
elseif strcmp(shape, 'parabolic_cylinder')
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = -2*x1 + x2 + x1^2 - 2*x1*x2 + x2^2;
    end
    
elseif strcmp(shape, 'wolverine')
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = 0.1 * (2*x1^2/x2 + exp(x2));
    end
    
elseif strcmp(shape, 'exponential')
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = exp(x1) + x2;
    end

elseif strcmp(shape, 'chair')
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = 0.5e-3 * (exp(x1) + x2^4);
    end
    
else
    error('Shape should be either "trough" or "paraboloid".')
end

end

