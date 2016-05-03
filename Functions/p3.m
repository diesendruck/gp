function [x, y] = p3(x, y, varargin)
% Shortcut to scatterplot, for when you have a matrix of x, and response y.
%
% Args:
%   x: n x d matrix of data values.
%   y: n x 1 matrix of response values.
%   mesh_gran: Number of ticks on mesh for plotting.
%
% Returns:
%   x: n x d matrix of data values.
%   y: n x 1 matrix of response values.

%scatter3(x(:,1), x(:,2), y);

if nargin > 2
    mesh_gran = varargin{1};
else
    mesh_gran = 2*sqrt(length(x));
end

%figure; 
do_buffer = 0;
[~, ~, ~, ~, ~, ~, xt1, xt2, xt] = compute_mesh_info(x, mesh_gran, ...
    do_buffer);
yq = griddata(x(:, 1), x(:, 2), y, xt1, xt2);
mesh(xt1, xt2, yq);
hold on;
scatter3(x(:,1), x(:,2), y);


end

