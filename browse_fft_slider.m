function [] = browse_fft_slider(X,Y,V,LS,I)
S = SpectroscopyData();
S.hf1 = figure;
S.hf1.Position = [2363 292 1146 530];
S.p = panel();
S.p.pack(1,2);
S.p.de.margin = 20;
S.p.margin = [15 25 10 10]; %left, bottom, right, top
% S.p.select('all');
% S.p.identitfy();

S.X = X;
S.Y = Y;
S.V = V;
S.LS = LS;
S.I = I;
S.X_cropped = X;
S.Y_cropped = Y;
S.LS_cropped = LS;
S.ctr = 0;

[S.hi1, S.ax1] = plot_real_space(S,1);
[S.hi2, S.ax2] = plot_fourier_space(S,1);

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
menu_radial_fft = uimenu(menu_main, 'Label', 'Radial FFT');
menu_radial_fft.Callback = @(hObject, eventData) radial_fft_call(hObject, eventData, S);
menu_make_video = uimenu(menu_main, 'Label', 'Make Video');
menu_make_video.Callback = @(hObject, eventData) makevideo(hObject, eventData, S);

% Update custom datatip
dcm_obj = datacursormode(S.hf1);
dcm_obj.UpdateFcn = @(hObject, eventData) dcm_update_fcn(hObject, eventData, size(S.hi1.CData), [S.X(1) S.X(end)], [S.Y(1) S.Y(end)]);
end

function [hi1, ax1] = plot_real_space(S,ii)
% Plots energy slice in real space
ax1 = S.p(1,1).select();
ax1.Box = 'on';
LS_realspace_slice = squeeze(S.LS(:,:,ii));
hi1 = imagesc(S.X,S.Y,LS_realspace_slice);
% axis([S.X(1) S.X(end) S.Y(1) S.Y(end)]);
axis tight;
set(gca,'YDir','normal');
colormap(ax1, parula);
colorbar();
[cmin, cmax] = color_scale(LS_realspace_slice, 3);
caxis(ax1, [cmin cmax]);
xlabel('X (nm)','FontSize',14);
ylabel('Y (nm)','FontSize',14);
title(ax1, ['E = ' num2str(S.V(ii)) ' eV'], 'fontsize', 14);
end

function LS_fft = calculate_fourier(LS)
% Computes the Fourier transform
LS_fft = fftshift(fft2(LS, 1.*size(LS,1), 1.*size(LS,2)));
LS_fft = abs(LS_fft);
% LS_fft = log(LS_fft);
% LS_fft = angle(LS_fft);
% LS_fft = imgaussfilt(LS_fft, 1);

% Remove dc by interpolation
dc_index = ceil((size(LS_fft,2)+1)/2);
LS_fft_without_dc = LS_fft;
LS_fft_without_dc(:,dc_index) = [];
% interp_dc_value = interp2(LS_fft_without_dc,(dc_index-1)/2,1:size(LS_fft,1), 'nearest');
% LS_fft(:,dc_index) = interp_dc_value;
LS_fft = LS_fft_without_dc;

% Symmeterization
% LS_fft = (LS_fft + imrotate(LS_fft, 120, 'crop') + imrotate(LS_fft, 240, 'crop'))/3;
end

function [hi2, ax2] = plot_fourier_space(S,ii)
% Computes the Fourier transform 
LS_fft = calculate_fourier(squeeze(S.LS(:,:,ii)));

LX = length(S.X);
dX = mean(diff(S.X));
q_x = pi.*linspace(-1,1,LX).*(1/dX);
LY = length(S.Y);
dY = mean(diff(S.Y));
q_y = pi.*linspace(-1,1,LY).*(1/dY);

% Plots the FT
ax2 = S.p(1,2).select();
ax2.Box = 'on';
hi2 = imagesc(q_x, q_y, LS_fft);
% axis([q_x(1) q_x(end) q_y(1) q_y(end)]);
axis image;
set(gca,'YDir','normal');
colormap(ax2, parula);
colorbar();
[cmin, cmax] = color_scale(LS_fft, 3);
caxis(ax2, [cmin cmax]);
xlabel('q_x (nm^{-1})','FontSize',12);
ylabel('q_y (nm^{-1})','FontSize',12);
title(ax2, ['E = ' num2str(S.V(ii)) ' eV'], 'fontsize', 14);
end
%% Create Functions
function [] = slider_CreateFcn(hObject,eventData,S)
addlistener(hObject, 'ContinuousValueChange', @(hObject,eventData)slider_call(hObject,eventData,S));
end

%% Callback functions
function [] = slider_call(varargin)
% Callback for the slider.
[h,S] = varargin{[1,3]};  % calling handle and data structure.
LS_realspace = squeeze(S.LS(:,:,round(get(h,'value'))));
I_slice = squeeze(S.I(:,:,round(get(h,'value'))));
% LS_realspace = LS_realspace./I_slice;
% LS_realspace = bsxfun(@minus, LS_realspace, mean(LS_realspace,2));
% LS_realspace = imgaussfilt(LS_realspace,1);
set(S.hi1,'cdata', LS_realspace);
[cmin, cmax] = color_scale(LS_realspace, 3);
caxis(S.ax1, [cmin cmax]);
title(S.ax1, ['E = ', sprintf('%0.3d',round(S.V(round(get(h,'value')))*1e3)), ' meV'], 'fontsize', 14);

LS_fft = calculate_fourier(squeeze(S.LS_cropped(:,:,round(get(h,'value')))));
set(S.hi2,'cdata',LS_fft);
[~, cmax] = color_scale(LS_fft, 3);
caxis(S.ax2, [0 cmax]);
title(S.ax2, ['E = ', sprintf('%0.3d',round(S.V(round(get(h,'value')))*1e3)), ' meV'], 'fontsize', 14)
drawnow;

assignin('base', 'S', S);

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
crop_vector = round(pos);
x = find(S.X>=crop_vector(1) & S.X<=crop_vector(1)+crop_vector(3));
y = find(S.Y>=crop_vector(2) & S.Y<=crop_vector(2)+crop_vector(4));
S.X_cropped = S.X(x);
S.Y_cropped = S.Y(y);
S.LS_cropped = S.LS(y,x,:);

slider_call(S.slider,0,S);
end

function [] = reset_button_call(varargin)
% Callback for the reset button.
S = varargin{3};  % calling handle and data structure.
S.X_cropped = S.X;
S.Y_cropped = S.Y;
S.LS_cropped = S.LS;
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
v = VideoWriter('MapMovie.avi');
v.FrameRate = 1;
open(v);
for i = 1:numel(S.V)
    S.slider.Value = i;
    slider_call(S.slider, [], S);
%     export_fig(S.hf1, sprintf('img/%0.3d.png', round(i)), '-png', '-nocrop', '-transparent', '-q90', '-r100');
    img = export_fig(S.hf1, '-png', '-nocrop', '-transparent', '-q90', '-r100');
    writeVideo(v, img);
end
close(v);
end