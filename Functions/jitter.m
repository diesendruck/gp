function [x_nsy_jit] = jitter(x_nsy)
% Jitter points, ensuring identical boundary after jitter.
%
% Args:
%   x_nsy: Matrix of convex+noise data points.
%
% Returns:
%   x_nsy_jit: Jittered inputs.

x_nsy_jit = x_nsy;

x1_l = min(x_nsy(:, 1)); x1_h = max(x_nsy(:, 1));
x2_l = min(x_nsy(:, 2)); x2_h = max(x_nsy(:, 2));
x1_range = x1_h - x1_l; x2_range = x2_h - x2_l;

% Jitter all but first and last points, which establish the boundary.
for i = 2:(length(x_nsy)-1)
    x1 = x_nsy(i, 1); x2 = x_nsy(i, 2);
    
    % For each component, if min, jitter up; if max, jitter down; otherwise
    % jitter randomly.
    if x1 == x1_l
        x_nsy_jit(i, 1) = x_nsy(i, 1) + abs((1e-6)*x1_range*randn());
    elseif x1 == x1_h
        x_nsy_jit(i, 1) = x_nsy(i, 1) - abs((1e-6)*x1_range*randn());
    else
        x_nsy_jit(i, 1) = x_nsy(i, 1) + (1e-6)*x1_range*randn();
    end
    
    if x2 == x2_l
        x_nsy_jit(i, 2) = x_nsy(i, 2) + abs((1e-6)*x2_range*randn());
    elseif x2 == x2_h
        x_nsy_jit(i, 2) = x_nsy(i, 2) - abs((1e-6)*x2_range*randn());
    else
        x_nsy_jit(i, 2) = x_nsy(i, 2) + (1e-6)*x2_range*randn();
    end
    
end

end

