function [cmin, cmax] = color_scale(data, sigma)
%COLOR_SCALE determines the color scale of data
%   COLOR_SCALE(DATA, SIGMA) Data is MxN 2D matrix. Returns cmin
%   and cmax.
%   Copyright 2016 WASP. 

% Overview of the algorithm:
%    Finds the mean and std of the 2d data.
%    cmin = mean-sigma*std;
%    cmax = mean+sigma*std;

if ~exist('data', 'var')
    disp('I am not a genie');
    return;
end
if ~exist('sigma', 'var')
    sigma = 3;
end

data_mean = mean(data(:), 'omitnan');
data_std = std(data(:), 'omitnan');
cmin = data_mean-sigma*data_std;
cmax = data_mean+sigma*data_std;
end