function [ytruth_on_mcmcgrid] = compute_truth_from_xt(xt, shape)
% Evaluate true convex function on xt grid (from run_gpmc).
% NOTE1: Functions appear in BOTH this file and make_noisy_convex.m.
% NOTE2: Each function gets a unique error variance, sig.
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
        ytruth_on_mcmcgrid(ii) = 0.2*(x1-3)^2;            
    end
    
elseif strcmp(shape, 'paraboloid')
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = 0.2 * ((x1-5)^2 + (x2-5)^2);
    end
    
elseif strcmp(shape, 'hand')
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = 1e-5*x1^6 - log(100*x2+1000) + 8;
    end
    
elseif strcmp(shape, 'parabolic_cylinder')
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = 0.1*(-2*x1 + x2 + x1^2 - 2*x1*x2 + x2^2);
    end
    
elseif strcmp(shape, 'wolverine')
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = 0.1 * (2*(x1-5)^2/(x2+1) + exp(2*(x2+1)/5));
    end
    
elseif strcmp(shape, 'exponential')
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = 0.3*(exp(x1/3) + 0.3*x2);
    end

elseif strcmp(shape, 'chair')
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = 1e-3 * (exp(x1/1.1) + 2*(x2-5)^4);
    end

elseif strcmp(shape, 'hannah2')
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = (x1 + x2)^2;
    end

elseif strcmp(shape, 'cm1')
    for ii = 1:n
        x1 = xt(ii, 1); x2 = xt(ii, 2);
        ytruth_on_mcmcgrid(ii) = 0.025 * (x1 + x2)^2;
    end
    
else
    error('Shape not recognized.')
end

end

