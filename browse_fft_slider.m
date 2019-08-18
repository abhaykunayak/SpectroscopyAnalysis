function [] = browse_fft_slider(X,Y,V,LS,I,Z)
S = SpectroscopyData();
S.hf1 = figure;
S.hf1.Position = [2363 292 1146 530];

S.X = X;
S.Y = Y;
S.Z = Z;
S.V = V;
S.LS = LS;
S.LS_fft = [];
S.I = I;
S.crop_pos = [];
S.I_cropped = I;
S.X_cropped = X;
S.Y_cropped = Y;
S.LS_cropped = LS;
S.Z_cropped = Z;
S.ctr = 0;

try
    S = plot_real_space(S);
    S = plot_fourier_space(S);
catch
end

S.slider = uicontrol('style','slide',...
                 'min',1,'max',numel(S.V),'val',1,...
                 'position',[400 10 300 30],...
                 'sliderstep',[1/numel(S.V) 1/numel(S.V)],...
                 'createfcn', {@slider_CreateFcn,S});
S.crop_button = uicontrol('Style', 'pushbutton', 'String', 'CROP',...
        'Position', [200 10 50 30],...
        'Callback', {@crop_button_call,S});
S.reset_button = uicontrol('Style', 'pushbutton', 'String', 'RESET',...
        'Position', [270 10 50 30],...
        'Callback', {@reset_button_call,S});

% Menu - Analysis
menu_main = uimenu(S.hf1, 'Label', 'Analysis');
menu_line_profile = uimenu(menu_main, 'Label', 'Line Profile');
menu_line_profile.Callback = @(hObject, eventData) line_profile_call(hObject, eventData, S);
menu_energy_profile = uimenu(menu_main, 'Label', 'Energy Profile');
menu_energy_profile.Callback = @(hObject, eventData) energy_profile(hObject, eventData, S);
menu_radial_fft = uimenu(menu_main, 'Label', 'Radial FFT');
menu_radial_fft.Callback = @(hObject, eventData) radial_fft_call(hObject, eventData, S);
menu_make_video = uimenu(menu_main, 'Label', 'Make Video');
menu_make_video.Callback = @(hObject, eventData) makevideo(hObject, eventData, S);

% Update custom datatip
dcm_obj = datacursormode(S.hf1);
dcm_obj.UpdateFcn = @(hObject, eventData) dcm_update_fcn(hObject, eventData, size(S.hi1.CData), [S.X(1) S.X(end)], [S.Y(1) S.Y(end)]);
end

function S = plot_real_space(S)
% Plots energy slice in real space
S.ax1 = subplot(1,2,1);
S.ax1.Box = 'on';
S = calculate_realspace(S, 1);
S.hi1 = imagesc(S.X,S.Y,S.LS_realspace);
axis image;
set(S.ax1, 'YDir', 'normal');
set(S.ax1, 'Layer', 'Top');
colormap(S.ax1, parula);
colorbar();
[cmin, cmax] = color_scale(S.LS_realspace, 3);
caxis(S.ax1, [cmin cmax]);
xlabel('X (nm)','FontSize',14);
ylabel('Y (nm)','FontSize',14);
title(S.ax1, ['E = ' num2str(S.V(1)) ' eV'], 'fontsize', 14);
end

function S = plot_fourier_space(S)
% Computes the Fourier transform 
S = calculate_fourier(S);

LX = length(S.X);
dX = mean(diff(S.X));
S.qx = pi.*linspace(-1,1,LX).*(1/dX);
LY = length(S.Y);
dY = mean(diff(S.Y));
S.qy = pi.*linspace(-1,1,LY).*(1/dY);

% Plots the FT
S.ax2 = subplot(1,2,2);
S.ax2.Box = 'on';
S.hi2 = imagesc(S.qx, S.qy, abs(S.LS_fft));
axis image;
set(S.ax2, 'YDir', 'normal');
set(S.ax2, 'Layer', 'Top');
colormap(S.ax2, flipud(gray));
colorbar();
[cmin, cmax] = color_scale(abs(S.LS_fft), 3);
caxis(S.ax2, [cmin cmax]);
xlabel('q_x (nm^{-1})','FontSize',12);
ylabel('q_y (nm^{-1})','FontSize',12);
title(S.ax2, ['E = ' num2str(S.V(1)) ' eV'], 'fontsize', 14);
end
%% Create Functions
function [] = slider_CreateFcn(hObject,eventData,S)
addlistener(hObject, 'ContinuousValueChange', @(hObject,eventData)slider_call(hObject,eventData,S));
end

%% Callback functions
function [] = slider_call(varargin)
% Callback for the slider.
[h,S] = varargin{[1,3]};  % calling handle and data structure.

% Update realspace spectroscopy
S = update_realspace(h,S);

% Update cropper
himr = findobj('Tag', 'imrect');
if isgraphics(himr)
    crop_data(S.crop_pos, S);
else
    S.LS_realspace_cropped = S.LS_realspace;
end

% Update fourier space spectroscopy
S = update_fourier_space(h,S);

% Line profile update
% calls the update function only if the line profile analysis figure is
% open
if isfield(S.data, 'ax') && ishandle(S.data.ax)
    line_profile_run(round(get(h,'value')),0,S);
end

% Radial FFT update
% calls the update function only if the radial fft figures are
% open
if isfield(S.data, 'hi_unfolded_qpi') && ishandle(S.data.hi_unfolded_qpi)
    update_radial_fft(round(get(h,'value')),0,S);
end

assignin('base', 'S', S);
end

function S = calculate_realspace(S, m)
% Manipulate the realspace
S.LS_realspace = squeeze(mean(S.LS(:,:,m),3));
S.I_slice = squeeze(mean(S.I(:,:,m),3));
% S.LS_realspace = tanh(S.LS_realspace./S.I_slice);
% S.LS_realspace = bsxfun(@minus, S.LS_realspace, smooth(mean(S.LS_realspace,2),50));
% S.LS_realspace = imrotate(S.LS_realspace, 33, 'bicubic', 'crop');
% S.LS_realspace = diff(S.LS_realspace,1,1);
% [FX, FY] = gradient(S.LS_realspace);
% S.LS_realspace = sqrt(FX.^2+FY.^2);
% S.LS_realspace = imgaussfilt(S.LS_realspace, 0.5);
end

function S = calculate_fourier(S)
% Computes the Fourier transform
% LS_slice = imgaussfilt(S.LS_realspace, 10);
LS_slice = S.LS_realspace_cropped;
LS_slice_sz = size(LS_slice);
wt1 = tukeywin(LS_slice_sz(1),0.3);
wt2 = tukeywin(LS_slice_sz(2),0.3);
LS_slice = LS_slice.*(wt1*wt2');
S.LS_fft = fftshift(fft2(LS_slice,...
    1.*LS_slice_sz(1), 1.*LS_slice_sz(2)));
% S.LS_fft = angle(S.LS_fft);
% S.LS_fft = unwrap(angle(S.LS_fft),2*pi);
S.LS_fft = abs(S.LS_fft);


% Z fft
% wt1 = tukeywin(size(S.Z,1),0.3);
% wt2 = tukeywin(size(S.Z,2),0.3);
% Z = S.Z.*(wt1*wt2');
% Z = S.Z;
% Z_fft = abs(fftshift(fft2(Z)));
% Z_fft = imgaussfilt(Z_fft, 0.5);
% S.LS_fft = S.LS_fft./Z_fft;

% S.LS_fft = log2(S.LS_fft);
S.LS_fft = imgaussfilt(S.LS_fft, 0.5);

% Remove dc by interpolation
% dc_index = ceil((size(LS_fft,2)+1)/2);
% LS_fft_without_dc = LS_fft;
% LS_fft_without_dc(:,dc_index) = [];
% interp_dc_value = interp2(LS_fft_without_dc,(dc_index-1)/2,1:size(LS_fft,1), 'nearest');
% LS_fft(:,dc_index) = interp_dc_value;
% LS_fft = LS_fft_without_dc;

% Symmeterization
% T = [1.2 0.08 0;
%      0 1 0;
%      0 0 1];
% tform = affine2d(T);
% LS_fft = imwarp(LS_fft, tform);
% LS_fft = 0*LS_fft + 1*imrotate(LS_fft, 117, 'crop');
% LS_fft = (LS_fft + imrotate(LS_fft, 120, 'crop') + imrotate(LS_fft, 240, 'crop'))/3;

end

function S = update_realspace(h,S)
% Updates the realspace.
n = round(get(h,'value'));
% m = n:(n+5);
m = n;
S = calculate_realspace(S, m);
set(S.hi1,'cdata', S.LS_realspace);
[cmin, cmax] = color_scale(S.LS_realspace, 3);
caxis(S.ax1, [cmin cmax]);
title(S.ax1, ['E = ', sprintf('%0.3f',S.V(round(get(h,'value')))*1e3), ' meV'], 'fontsize', 14);
end

function S = update_fourier_space(h,S)
% Updates the fourier space.
% S.LS_fft = calculate_fourier(squeeze(S.LS_cropped(:,:,round(get(h,'value')))));
% LS_realspace = squeeze(S.LS_cropped(:,:,round(get(h,'value'))));
n = round(get(h,'value'));
% S = calculate_realspace(S, n);
S = calculate_fourier(S);

LX = length(S.X_cropped);
dX = mean(diff(S.X_cropped));
S.qx = pi.*linspace(-1,1,LX).*(1/dX);
LY = length(S.Y_cropped);
dY = mean(diff(S.Y_cropped));
S.qy = pi.*linspace(-1,1,LY).*(1/dY);

set(S.hi2, 'xdata', S.qx);
set(S.hi2, 'ydata', S.qy);
% set(S.hi2, 'cdata', abs(S.LS_fft));
set(S.hi2, 'cdata', S.LS_fft);

[cmin, cmax] = color_scale(abs(S.LS_fft), 2);
% caxis(S.ax2, [cmax*0.2 cmax*1.5]);
% caxis(S.ax2, [-pi pi]);
caxis(S.ax2, 'auto');
title(S.ax2, ['E = ', sprintf('%0.3f',S.V(round(get(h,'value')))*1e3), ' meV'], 'fontsize', 14);
end

function [] = crop_button_call(varargin)
% Callback for the crop button.
S = varargin{3};  % calling handle and data structure.

hi = imrect(S.ax1);
func_new_pos = @(p) crop_data(p,S);
addNewPositionCallback(hi,func_new_pos);
fcn = makeConstrainToRectFcn('imrect',get(S.ax1,'XLim'),get(S.ax1,'YLim'));
setPositionConstraintFcn(hi,fcn);
drawnow;
end

function [] = crop_data(pos,S)
% Crops data
S.crop_pos = pos;
crop_vector = round(pos);
x = find(S.X>=crop_vector(1) & S.X<=crop_vector(1)+crop_vector(3));
y = find(S.Y>=crop_vector(2) & S.Y<=crop_vector(2)+crop_vector(4));
S.X_cropped = S.X(x);
S.Y_cropped = S.Y(y);
S.LS_cropped = S.LS(y,x,:);
S.I_cropped = S.I(y,x,:);
S.LS_realspace_cropped = S.LS_realspace(y,x);
S.I_slice_cropped = S.I_slice(y,x);
% S.Z_cropped = S.Z(y,x);

% slider_call(S.slider,0,S);
update_fourier_space(S.slider,S);
end

function [] = reset_button_call(varargin)
% Callback for the reset button.
S = varargin{3};  % calling handle and data structure.
S.X_cropped = S.X;
S.Y_cropped = S.Y;
S.LS_cropped = S.LS;
S.I_cropped = S.I;
S.LS_realspace_cropped = S.LS_realspace;
S.I_slice_cropped = S.I_slice;
delete(findobj('Tag', 'imrect'));
slider_call(S.slider,0,S);
end

function [] = line_profile_call(varargin)
% Callback for menu: Line Profile
% Setups the GUI for analysis
S = varargin{3};

fig = figure;
fig.Position = [1134 179 560 400];
ax = axes;
ax.Box = 'on';
ax.Position = [0.1300 0.3200 0.7750 0.6000];
ax.Title.String = 'Line Profile Analysis';
ax.YLabel.String = 'dI/dV (au)';
ax.XLabel.String = 'Position';

panel_main = uipanel(fig);
panel_main.Position = [0 0 1.0000 0.2075];
panel_main.Title = 'Options';

txtbox_x = uicontrol(panel_main,'Style','text',...
                'String','Center X',...
                'Position',[20 40 60 20]);
inputbox_x = uicontrol(panel_main,'Style','edit',...
                'String',num2str(size(S.hi1.CData,1)/2),...
                'Position',[80 40 60 20]);            
txtbox_y = uicontrol(panel_main,'Style','text',...
                'String','Center Y',...
                'Position',[140 40 60 20]);
inputbox_y = uicontrol(panel_main,'Style','edit',...
                'String',num2str(size(S.hi1.CData,2)/2),...
                'Position',[200 40 60 20]);            
txtbox_theta = uicontrol(panel_main,'Style','text',...
                'String','theta',...
                'Position',[20 20 60 20]);
inputbox_theta = uicontrol(panel_main,'Style','edit',...
                'String','0',...
                'Position',[80 20 60 20]);            
txtbox_radius = uicontrol(panel_main,'Style','text',...
                'String','radius',...
                'Position',[140 20 60 20]);
inputbox_radius = uicontrol(panel_main,'Style','edit',...
                'String','50',...
                'Position',[200 20 60 20]);
txtbox_delta_theta = uicontrol(panel_main,'Style','text',...
                'String','d theta',...
                'Position',[260 20 60 20]);
inputbox_delta_theta = uicontrol(panel_main,'Style','edit',...
                'String','15',...
                'Position',[320 20 60 20]);
btn_run = uicontrol(panel_main, 'Style', 'pushbutton', ...
                'String', 'Run',...
                'Position', [440 40 50 20],...
                'Callback', {@line_profile_run,S}); 
btn_clear = uicontrol(panel_main, 'Style', 'pushbutton', ...
                'String', 'Clear',...
                'Position', [440 20 50 20],...
                'Callback', 'cla');

% Add the relevan fields to S.data which is to accesed in other parts of
% the program
S.data.ax = ax;
S.data.inputbox_x = inputbox_x;
S.data.inputbox_y = inputbox_y;
S.data.inputbox_theta = inputbox_theta;
S.data.inputbox_radius = inputbox_radius;
S.data.inputbox_delta_theta = inputbox_delta_theta;
end

function [] = line_profile_run(varargin)

% Check if function is being called by the slider, if yes then set the
% offset to non zero values. This will help to better visualize the plots
% for different energies. 
if isnumeric(varargin{1})
    offset = varargin{1}-1;
else
    offset = 0;
end
S = varargin{3};

% Get the relevant data from the input text boxes
center_x = str2double(S.data.inputbox_x.String);
center_y = str2double(S.data.inputbox_y.String);
start_x = center_x;
start_y = center_y;

theta = str2double(S.data.inputbox_theta.String);
radius = str2double(S.data.inputbox_radius.String);
avg_over_theta = str2double(S.data.inputbox_delta_theta.String);
LS_E = S.hi1.CData;

% Line profile of single pixel wide line
% end_x = center_x + round(radius*cosd(theta));
% end_y = center_y + round(radius*sind(theta));
% LS_E = S.hi1.CData;
% xi = [start_x end_x];
% yi = [start_y end_y];
% [cx, cy, c] = improfile(LS_E, xi, yi, 'nearest');

% Line profile of multi pixel wide line, all same angle
% avg_over_lines = 5;
% for i = 1:avg_over_lines
%     if i>(avg_over_lines/2)
%         j = round((avg_over_lines/2))-i;
%     else
%         j = i;
%     end
%     start_x = center_x - round((j-1)*sind(theta));
%     start_y = center_y + round((j-1)*cosd(theta));
%     end_x = start_x + round(radius*cosd(theta));
%     end_y = start_y + round(radius*sind(theta));
%     xi(:,i) = [start_x end_x];
%     yi(:,i) = [start_y end_y];
%     c(:,i) = improfile(LS_E, xi(:,i), yi(:,i), 'nearest');
% end

% Line profile averaged over an angle
[xi, yi, c] = line_profile(theta, avg_over_theta, start_x, start_y, radius, LS_E, offset);
avg_c = mean(squeeze(c),2);

% Plot the line profile
hold(S.data.ax, 'on');
plot(S.data.ax, avg_c, 'LineWidth', 2);
hold(S.data.ax, 'off');

% Plot lines to indicate the range over which line profiles have been
% averaged over
hold(S.ax1, 'on');
line(S.ax1, S.X(xi(:,1)), S.Y(yi(:,1)), 'LineWidth', 1, 'Color', [1 1 1]);
line(S.ax1, S.X(xi(:,avg_over_theta-1)), S.Y(yi(:,avg_over_theta-1)), ...
    'LineWidth', 1, 'LineStyle', '--', 'Color', [1 0 0]);
line(S.ax1, S.X(xi(:,end)), S.Y(yi(:,end)), 'LineWidth', 1, ...
    'LineStyle', '--', 'Color', [1 0 0]);
hold(S.ax1, 'off');

end

function [xi, yi, c] = line_profile(theta_vec, avg_over_theta, start_x, start_y, radius, LS_E, offset)
    for k = 1:numel(theta_vec)
        theta = theta_vec(k);
        for i = 1:2*avg_over_theta-1
            if i>avg_over_theta
                j = avg_over_theta-i;
            else
                j = i;
            end
            theta_i = theta + (j-1);
            end_x = start_x + round(radius*cosd(theta_i));
            end_y = start_y + round(radius*sind(theta_i));
            xi(:,i) = [start_x end_x];
            yi(:,i) = [start_y end_y];
            c(:,i,k) = improfile(LS_E+(0.5*offset), xi(:,i), yi(:,i), radius, 'nearest');
        end
    end
end

function [] = radial_fft_call(varargin)
S = varargin{3};
% Get the relevant data from the input text boxes
% center_x = 65;
% center_y = 78;
center_x = 71;
center_y = 87;
start_x = center_x;
start_y = center_y;

avg_over_theta = 5;
theta = 0:1:360;
radius = 64;
LS_E = S.hi1.CData;
offset = 0;

[~, ~, c] = line_profile(theta, avg_over_theta, start_x, start_y, radius, LS_E, offset);
c_real_space = squeeze(mean(c,2));
x = S.X(start_x:(start_x+radius));
x = x - min(x);

figure;
axes;
S.data.hi_unfolded_qpi = imagesc(x,theta,c_real_space');
set(gca, 'YDir', 'Normal');
xlabel('Position');
ylabel('\theta in deg');
title('Unfolded QPI');
caxis([-0.5 0.5]);

c_fft = abs(fftshift(fft(c_real_space, 1.*size(c_real_space,1), 1)));
LX = numel(x);
dX = mean(diff(x));
q_x = pi.*linspace(-1,1,LX).*(1/dX);

figure;
axes;
S.data.hi_unfolded_fft = imagesc(q_x, theta, c_fft');
set(gca, 'YDir', 'Normal');
xlabel('q');
ylabel('\theta in deg');
title('Unfolded FFT');
caxis([0 5]);
end

function [] = update_radial_fft(varargin)
S = varargin{3};
center_x = 71;
center_y = 87;
start_x = center_x;
start_y = center_y;

avg_over_theta = 5;
theta = 0:1:360;
radius = 64;
LS_E = S.hi1.CData;
offset = 0;

[~, ~, c] = line_profile(theta, avg_over_theta, start_x, start_y, radius, LS_E, offset);
c_real_space = squeeze(mean(c,2));
S.data.hi_unfolded_qpi.CData = c_real_space';

c_fft = abs(fftshift(fft(c_real_space, 1.*size(c_real_space,1), 1)));
S.data.hi_unfolded_fft.CData = c_fft';
end

function output_txt = dcm_update_fcn(obj,event_obj,size_data,x_cord,y_cord)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
output_txt = {['X: ',num2str(pos(1),4), ', ',num2str(axes2pix(size_data(1), x_cord, pos(1)))],...
    ['Y: ',num2str(pos(2),4), ', ',num2str(axes2pix(size_data(2), y_cord, pos(2)))]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end
end

function [] = makevideo(varargin)
S = varargin{3};
% v = VideoWriter('MapMovie.avi');
% v.FrameRate = 1;
% open(v);
for i = 1:numel(S.V)
    S.slider.Value = i;
    slider_call(S.slider, [], S);
    export_fig(S.hf1, sprintf('img/%0.3d.png', round(i)),...
        '-png', '-nocrop', '-transparent', '-q90', '-r100');
%     img = export_fig(S.hf1, '-png', '-nocrop', '-transparent', '-q90', '-r100');
%     writeVideo(v, img);
end
% close(v);
end

function [] = energy_profile(varargin)
% Callback for menu: Energy Profile
S = varargin{3};
figure;
ax = axes;
% theta = 0;
for theta=0:10:180
LS_lc = [];

for i = 1:numel(S.V)
    S.slider.Value = i;
    slider_call(S.slider, [], S);
    LS_rot = imrotate(S.LS_fft, theta, 'bicubic', 'crop');
    sz_lc = size(LS_rot);
    xidx = round(sz_lc(1)/2)+[-1:1];
    LS_lc(i,:) = squeeze(mean(LS_rot(xidx,:),1));
end

% Plot the energy profile

imagesc(ax, S.qy, S.V*1000, LS_lc);
% imagesc(S.qy, S.V*1000, LS_lc.*sign(S.qy));
set(ax, 'YDir', 'normal');
set(ax, 'Layer', 'Top');
colormap(ax, flipud(gray));
colorbar();
[~, cmax] = color_scale(LS_lc, 2);
caxis(ax, [cmax*0.1 cmax*1]);
% caxis(ax, [0 20]);
xlabel('q_x (nm^{-1})','FontSize',12);
ylabel('Energy (meV)','FontSize',12);
title(ax, ['theta = ', num2str(theta)], 'fontsize', 14);

% Save data to workspace
S.data = LS_lc;
pause;
end

% for i = 1:numel(S.V)
%     S.slider.Value = i;
%     slider_call(S.slider, [], S);
%     cdw_pos = [46,173; 73,207; 103,184; 105,128; 79,96; 49,118; 
%         58,197; 91,201;];
%     for j=1:size(cdw_pos,1)
%         p(j,i) = mean(mean(S.LS_fft((cdw_pos(j,1)-1):(cdw_pos(j,1)+1), ...
%             (cdw_pos(j,2)-1):cdw_pos(j,2)+1),1),2);
%     end
% 
% end
% 
% for j=1:size(cdw_pos,1)
%     hold on;
%     plot(ax, S.V, p(j,:), 'DisplayName', ['p',num2str(j)]);
%     hold off;
% end
% assignin('base', 'p', p);

end