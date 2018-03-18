% Reference 'LineSpec.m' by Jonathan. \MATLAB\STM Data Acquisition\LineSpec.m
% Read spectroscopy data from *.3ds file using Jonathan's modified code
% 'load3dsV.m'

%% Load .3ds file
[FILENAME, PATH] = uigetfile('*.3ds');
[X, Y, V, Z, DATA, header, ~] = load3dsV([PATH FILENAME]);
LS = squeeze(DATA(:,:,:,2));  % extracts the average channel of the line spectrum
I = squeeze(DATA(:,:,:,1));  % extracts the average current channel

%% Remove bad sweeps;
[LS, I] = remove_bad_sweeps(DATA, 2);

%% Fix data scale
LS = LS./max(LS(:));
I = I./max(I(:));
X = X*1e9;                      
Y = Y*1e9;
Z = Z*1e9;

% offs = 0;                   % Offset for clarity in nS
% sensitivity = 20e-3;        % LI7265 sensitivity
% scale = 0.1*sensitivity;
% dV = 3e-3;                  % OSC amplitude
% LS = LS*scale/dV;

%% Select lines
lines_selected = [1 2 3 4 5];
LS = LS(lines_selected,:,:);
I = I(lines_selected,:,:);
Z = Z(lines_selected,:);
Y = Y(lines_selected);

%% Subtract Linear fit
[Z, slope] = subtract_linear_fit(Z);                     

%% Correct Slope
slope_avg = mean(slope);
theta = atand(slope_avg);
X = X./cosd(theta);

%% Drift Correction
[Z, LS] = drift_correct(Z, LS);

%% Plot Topography 
plot_topography(X,Y,Z);

%% Data Analysis Each Energy SLice
browse_fft_slider(X,Y,V,LS,I);

%% Data Analysis Energy Dispersion
data_analysis(X,Y,Z,V,LS,I);
