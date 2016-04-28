function Pmono = monotone_2d(x_nsy, y_nsy)
% Project vector to monotone space.
%
% Args:
%   x_nsy: Matrix (nx2) of data input values.
%   y_nsy: Matrix (nx1) of noisy response values.
%
% Returns:
%   Pmono: Matrix (2x2) of the monotone response variable.

% Define convergence indicator and max number of iterations in while loop.
iter = 0;
max_iter = 1e3;
eps = 1e-6;

% Set up initial surface.
dim = sqrt(length(x_nsy));
x1 = reshape(x_nsy(:, 1), dim, dim)';
x2 = reshape(x_nsy(:, 2), dim, dim)';
y = reshape(y_nsy, dim, dim)';
f_init = y_nsy;

% First monotone projection, and its residual.
f_mono = monotone_1d(f_init);
r_mono = f_mono - f_init;

% First convex projection, and its residual.
f_conv = convex_1d(x_nsy, f_init + r_mono);
r_conv = f_conv - (f_init + r_mono);

while (iter < max_iter)
    
end
    
    
    
 

   


