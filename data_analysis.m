function [] = data_analysis(X,Y,Z,V,LS,I)
%DATA_ANALYSIS analysis tool for spectroscopic linecut data
%   DATA_ANALYSIS(X,Y,Z,V,LS) X,Y,Z are the spatial co-ordinates. V is the 
%   energy axis. LS is the 3-D spectroscopy data.
%   Copyright 2016 WASP. 

% Overview of the algorithm:
%    Plots topography, average point spectrum, avg real space QPI, 
%    FFT of real space averaged QPI and FFT of fourier space averaged QPI.

%% Define Structure
S = initialize_structure(X,Y,Z,V,LS,I);

%% Plot Topography
fig_topo = figure('Position', [550 560 600 400]);
% Dummy axes to just show X and Y in nm
topography_pos_axes = axes('Layer', 'bottom');
topography_pos_axes.Color = 'none';
axis([S.X(1) S.X(end) S.Y(1) S.Y(end)]);
xlabel('X (nm)','FontSize',14);
ylabel('Y (nm)','FontSize',14);
title('Topography','FontSize',12, 'VerticalAlignment', 'bottom');

% Axes which is used by imrect to do the cropping
topography_axes = axes;
imagesc(S.Z);
axis tight;
set(gca,'YDir','normal');
colormap(gray);
colorbar();
[cmin, cmax] = color_scale(S.Z, 2);
caxis([cmin cmax]);
topography_axes.XTickLabel = [];
topography_axes.YTickLabel = [];

topography_pos_axes.Position = topography_axes.Position;

%% Plot Point Spectroscopy
fig_ps = figure('Position', [1170 560 600 400]);
S.point_spectroscopy_axes = axes;
S.point_spectroscopy_plot = plot(S.V, S.LS_avg, 'LineWidth', 1.5);
xlabel('V (eV)');
ylabel('dI/dV (au)');
title('Average Point Spectroscopy');
axis tight;

%% Plot Energy vs X
fig_qpi = figure('Position', [550 80 600 400]);
S.energy_x_axes = axes;
S.energy_x_im = imagesc(S.X, S.V, S.LS_avg_map');
set(gca, 'YDir', 'Normal');
xlabel('X (nm)');
ylabel('E (eV)');
colormap();
colorbar;
title('QPI');
axis tight;
[cmin, cmax] = color_scale(S.LS_avg_map, 3);
caxis(S.energy_x_axes, [cmin cmax]);

%% Calculate FFT: averaged in real space
S = calc_fft(S);

%% Plot FFT: Energy vs q
fig_fft = figure('Position', [1170 80 600 400]);
S.energy_q_axes = axes;
S.energy_q_im = imagesc(S.q, S.V, S.LS_fft);
xlabel('q_{x} (nm^{-1})');
ylabel('E (eV)');
title('FFT');
axis tight;
set(gca,'YDir','normal');
colormap();
colorbar();
[~, cmax] = color_scale(S.LS_fft, 3);
caxis(S.energy_q_axes, [0 cmax]);

%% Align figure windows
iptwindowalign(fig_topo, 'right', fig_ps, 'left');
iptwindowalign(fig_topo, 'bottom', fig_qpi, 'top');
iptwindowalign(fig_topo, 'hcenter', fig_qpi, 'hcenter');
iptwindowalign(fig_ps, 'bottom', fig_fft, 'top');
iptwindowalign(fig_ps, 'hcenter', fig_fft, 'hcenter');

%% Select Image ROI
S.roi_imrect_handle = imrect(topography_axes);
func_new_pos = @(p) crop_data_roi(p,S);
addNewPositionCallback(S.roi_imrect_handle,func_new_pos);
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(S.roi_imrect_handle,fcn);

end

%% Initialize Structure
function S = initialize_structure(X,Y,Z,V,LS,I)
    S.X = X;
    S.X_cropped = X;
    S.Y = Y;
    S.Y_cropped = Y;
    S.Z = Z;
    S.V = V;
    S.LS = LS;
    S.LS_cropped = S.LS;
    S.I = I;
    S.x = 1:length(S.X);
    S.y = 1:length(S.Y);
    
    % Average Point Spectroscopy
    S.LS_avg = mean(squeeze(mean(S.LS,2)),1);
    S.LS_avg_cropped = S.LS_avg;
    
    % Average Energy vs X
    S.LS_avg_map = squeeze(mean(S.LS,1));
end

%% Calls calculate and update functions
function [] = crop_data_roi(pos,S)
% Crop Data
crop_vector = round(pos);
S.x = crop_vector(1):crop_vector(1)+crop_vector(3)-1;
S.y = crop_vector(2):crop_vector(2)+crop_vector(4)-1;
if length(S.x)<=length(S.X) && length(S.y)<=length(S.Y)
    S.X_cropped = S.X(S.x);
    S.Y_cropped = S.Y(S.y);
    S.LS_cropped = S.LS(S.y,S.x,:);
    S.I_cropped = S.I(S.y,S.x,:);
    S.Z_cropped = S.Z(S.y,S.x);
else
    warning('Crop dimension exceeds data dimension');
end

S = calc_avg_map(S);
S = calc_fft(S);

assignin('base', 'S', S);
update_data_roi(S);
end

%% Calculate Functions 
function S = calc_avg_map(S)
% Normalization
S = normalize_dIdV(S);

% Background subtraction
S = subtract_dIdV(S, 'subavg');

% Average in Y Energy vs X
S.LS_avg_map = squeeze(mean(S.LS_cropped,1));

% Derivative of average spectroscopy along Energy
% S.LS_avg_map = diff(S.LS_avg_map, 1, 1);

% Smooth average spectroscopy map
S.LS_avg_map = imgaussfilt(S.LS_avg_map, 0.5);
end

function S = normalize_dIdV(S)
% Normalization: Divide QPI by an energy integrated spectrum

% Integrate all energy
% S.LS_avg_map = bsxfun(@rdivide, S.LS_avg_map, mean(S.LS_avg_map,2));
for i=1:size(S.LS_cropped,1)
    S.LS_cropped(i,:,:) = bsxfun(@rdivide, squeeze(S.LS_cropped(i,:,:)), mean(squeeze(S.LS_cropped(i,:,:)),2));
end
% for i=1:size(S.LS_cropped,1)
% %     S.LS_cropped(i,:,:) = bsxfun(@rdivide, squeeze(S.LS_cropped(i,:,:)), smooth(squeeze(S.I_cropped(i,:,:))./S.V,10));
%     for j=1:size(S.LS_cropped,2)
%         S.LS_cropped(i,j,:) = squeeze(S.LS_cropped(i,j,:)).'./smooth(squeeze(S.I_cropped(i,j,:))./S.V.',1).';
% %         S.LS_cropped(i,j,:) = squeeze(S.LS_cropped(i,j,:)).'./abs(S.V);
%     end
% end

% Integrate energy from fermi energy to parking bias
% S.LS_avg_map = bsxfun(@rdivide, S.LS_avg_map, mean(S.LS_avg_map(:,65:128),2));
% for i=1:size(S.LS_cropped,1)
%     S.LS_cropped(i,:,:) = bsxfun(@rdivide, squeeze(S.LS_cropped(i,:,:)), mean(squeeze(S.LS_cropped(i,:,65:128)),2));
% end

% Integrate energy in a given range
% length_V = length(S.V);
% % norm_en_range = 10:20;
% norm_en_range = (length_V-15):(length_V-5);
% asymp_value = squeeze(mean(mean(S.LS_cropped(:,end-5:end,norm_en_range),2),3));
% no_norm_value = squeeze(mean(S.LS_cropped(:,:,norm_en_range),3));
% norm_factor = asymp_value./no_norm_value;
% for i=1:size(S.LS_cropped,2)
%     S.LS_cropped(1,i,:) = S.LS_cropped(1,i,:).*norm_factor(i);
% end
end

function S = subtract_dIdV(S, method)
% Background subtraction:

if strcmp(method, 'subavg')
    % Subtract average spectroscopy from the map
    for i=1:size(S.LS_cropped,1)
        S.LS_cropped(i,:,:) = bsxfun(@minus, squeeze(S.LS_cropped(i,:,:)), smooth(squeeze(mean(S.LS_cropped(i,:,:),2)),50)');
    end
elseif strcmp(method, 'subavg_reg')
    % Subtract average spectroscopy of a specific featureless region from the map
    S.LS_avg_map = bsxfun(@minus, S.LS_avg_map, squeeze(mean(S.LS(1,200:300,:),2))');
end

end

function S = calc_fft(S)
% Calculate FFT after Y averaged
LX = length(S.X_cropped);
dX = mean(diff(S.X_cropped));
S.q = pi.*linspace(-1,1,LX).*(1/dX);

% Fourier window: cos
% [M, N] = size(S.LS_avg_map);
% w = repmat(cos(linspace(-pi/2, pi/2, M)).', [1, N]);
% S.LS_avg_map = S.LS_avg_map.*w;

% Method 1: FT of the y averaged spectroscopy map 
% S.LS_fft = abs(fftshift(fft(S.LS_avg_map', LX,2),2));

% Method 2
LS_fft = abs(fftshift(fft(S.LS_cropped, LX,2),2));
S.LS_fft = squeeze(mean(LS_fft,1))';

% Fourier window: cos 2D
% [M, N, E] = size(S.LS_cropped);
% wm = cos(linspace(-pi/2, pi/2, M));
% wn = cos(linspace(-pi/2, pi/2, N));
% w = repmat(wm.' * wn, [1 1 E]);
% S.LS_cropped = S.LS_cropped.*w;

% Method 3: 2D FFT and then average over y
% LS_fft2 = abs(fftshift(fftshift(fft2(S.LS_cropped),2),1));
% S.LS_fft = squeeze(mean(LS_fft2,1)).';

% Remove dc peak
S.LS_fft = remove_dc(LX, S.LS_fft);

% Apply logscale
% S.LS_fft = log(S.LS_fft);
end

function LS_fft = remove_dc(LX, LS_fft)
% Replace dc value by average
dc_index = floor(LX/2)+1;
mean_index = [dc_index-2:dc_index-1, dc_index+1:dc_index+2];
LS_fft(:,dc_index) = mean(LS_fft(:,mean_index),2);
end

%% Calls all Update Functions
function [] = update_data_roi(S)
update_avg_spectroscopy(S);
update_avg_spectroscopy_map(S);
update_fft(S);
end

%% Update Funtions
function [] = update_avg_spectroscopy(S)
% Average Point Spectroscopy
S.LS_avg_cropped = squeeze(mean(mean(S.LS_cropped,2),1));

% Update Plot
set(S.point_spectroscopy_plot, 'XData', S.V);
set(S.point_spectroscopy_plot, 'YData', S.LS_avg_cropped);
title(S.point_spectroscopy_axes, ['Average Point Spectroscopy: Rows #', num2str(S.y(1)), '-', num2str(S.y(end))]);
end

function [] = update_avg_spectroscopy_map(S)
% Update Plot
set(S.energy_x_im, 'XData', S.X_cropped);
set(S.energy_x_im, 'YData', S.V);
set(S.energy_x_im, 'CData', S.LS_avg_map');
[cmin, cmax] = color_scale(S.LS_avg_map, 3);
caxis(S.energy_x_axes, [cmin/2 cmax/2]);
end

function [] = update_fft(S)
% Update Plot
set(S.energy_q_im, 'XData', S.q);
set(S.energy_q_im, 'YData', S.V);
set(S.energy_q_im, 'CData', S.LS_fft);
[~, cmax] = color_scale(S.LS_fft, 3);
caxis(S.energy_q_axes, [0 cmax]);
end