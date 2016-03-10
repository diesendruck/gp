function [x, y] = p3(x, y)
% Shortcut to scatterplot, for when you have a matrix of x, and response y.
%
% Args:
%   x: n x d matrix of data values.
%   y: n x 1 matrix of response values.
%
% Returns:
%   x: n x d matrix of data values.
%   y: n x 1 matrix of response values.
scatter3(x(:,1), x(:,2), y);

end

