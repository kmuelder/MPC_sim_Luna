% Read semicolon-delimited CSV file
filename = fullfile('input', 'Lab1_clean.csv');

% Read the data (assuming numeric values)
data = readmatrix(filename, 'Delimiter', ';');

% Remove empty columns if any
data(:, all(isnan(data), 1)) = [];

% Plot heatmap
figure;
imagesc(data);
colormap('hot'); % or 'jet', 'parula', etc.
colorbar;
axis equal tight;
title('Beam Profile Heatmap');
xlabel('X Pixels');
ylabel('Y Pixels');

% Enable interactive data cursor
datacursormode on;
