function [Z_subtracted, m] = subtract_linear_fit(Z)
%SUBTRACT_LINEAR_FIT subtracts the linear fit from Z
%   SUBTRACT_LINEAR_FIT(Z) LS is a 2d matrix data Z(X,Y). Returns the
%   subtracted Z data. Assumes X is uniform. For faster performance start
%   the parralel pool beforehand.
%   Copyright 2016 WASP. 

% Overview of the algorithm:
%    Performs linear fitting ('poly1') to each row of Z 
%    and then returns the subtracted row data. 

% Empty matrix for perfomance
X = 1:size(Z,2);
Z_subtracted = zeros(size(Z));
m = zeros(size(Z,1),1);
for i = 1:size(Z,1)
        % Prepare data to fit
        [xData, yData] = prepareCurveData(X, Z(i, :));

        % Set up fittype and options.
        ft = fittype('poly1');

        % Fit model to data.
        fitresult = fit(xData, yData, ft);

        % Subtract linear fit from data
        yData_fit = fitresult(xData);
        yData_subtracted = yData - yData_fit;

        % Create subtracted image
        Z_subtracted(i, :) = yData_subtracted;
        
        % Store slope fit parameter
        m(i) = fitresult.p1;
end
end