function [d] = get_shape_dimension(shape)
% Given shape name (String), return shape dimension d.

d = 0;
one_dim_shapes = {'parabola' 'flat_parabola'};
two_dim_shapes = {'trough' 'paraboloid' 'parabolic_cylinder'};

if strmatch(shape, one_dim_shapes)
    d = 1;
elseif strmatch(shape, two_dim_shapes)
    d = 2;
end

end

