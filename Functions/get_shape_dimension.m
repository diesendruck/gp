function [d] = get_shape_dimension(shape)
% Given shape name (String), return shape dimension d.
%
% Args:
%   shape: String indicating which convex function to use.
%
% Returns:
%   d: Integer number of dimensions of each data point.

d = 0;
one_dim_shapes = {'parabola' 'flat_parabola'};
two_dim_shapes = {'trough' 'paraboloid' 'parabolic_cylinder', 'wave', ...
    'plane', 'product', 'exponential'};

if strmatch(shape, one_dim_shapes)
    d = 1;
elseif strmatch(shape, two_dim_shapes)
    d = 2;
end

end

