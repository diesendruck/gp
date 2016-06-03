% Implement monotone on set of functions

% Import data.                           
filename = 'data/UnconstrainedFunctionEstimates_To_Mo_2016-06-01.csv';
xydat = importdata(filename);
x_nsy = xydat(:, 1);

y_samples = xydat(:, [2:size(xydat, 2)]);
num_samples = size(y_samples, 2);

results = zeros(size(xydat));
results(:, 1) = x_nsy;

% Compute each monotone function.
for i = 1:num_samples
    y_nsy = y_samples(:, i);
    
    y_mono = monotoneUp_1d(y_nsy);
    
    results(:, i+1) = y_mono;  % Plus one, because x_nsy is in first column.
end

% Show each result.
for j = 1:num_samples
    figure
    plot(y_samples(:, j), 'r.', 'markersize', 10); hold on;
    plot(results(:, j+1), 'b*'); hold on;
    plot(y_samples);
end

csvwrite('data/UnconstrainedFunctionEstimates_To_Mo_2016-06-01_MONOTONE.csv', results);
