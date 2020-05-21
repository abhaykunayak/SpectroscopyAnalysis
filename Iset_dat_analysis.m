%% Load multiple dat files
dirname = uigetdir;
datfiles = dir([dirname, '\*.dat']);
data = [];
info = struct;
Z = zeros(size(datfiles));
Iset = zeros(size(datfiles));
for i = 1:length(datfiles)
    disp([num2str(i), ' of ', num2str(length(datfiles))]);
    [info(i).header, data(:,i,:), info(i).channels] = ...
        loaddat(fullfile(dirname, datfiles(i).name));
    Z(i) = str2double(info(i).header.z_controller_z__m_);
    Iset(i) = str2double(info(i).header.z_controller_setpoint);
end

%% Select good sweeps
PS_all = data(:,:,6:3:end).*1e12;
I_all = data(:,:,5:3:end).*1e9;

PS_std = std(PS_all,[],3);
PS_mean = mean(PS_all,3);

PS_mask = double((PS_all-PS_mean)>(PS_std*0.5));
PS_mask(PS_mask==0) = nan;
PS = mean(PS_all.*PS_mask,3,'omitnan');
I = mean(I_all.*PS_mask,3,'omitnan');

V = squeeze(data(:,1,1)).*1e3;

Z = (Z(1)-Z).*1e9;
Iset = Iset.*1e9;

%%
PS_avg_fwd = squeeze(data(:,:,3)).*1e12;
I_avg_fwd = squeeze(data(:,:,2)).*1e9;
% PS_avg_bwd = squeeze(data(:,:,66)).*1e12;
% I_avg_bwd = squeeze(data(:,:,65)).*1e9;
V = squeeze(data(:,1,1)).*1e3;
% 
Z = (Z(1)-Z).*1e9;
Iset = Iset.*1e9;

%% Normaliza data
PS = PS_avg_fwd;
% PS = (PS_avg_fwd+PS_avg_bwd)/2;
% PS = PS_avg_fwd./mean(I_avg_fwd(1:2,:),1);
% PS = PS_avg_fwd./mean(PS_avg_fwd(1:10,:),1);
% PS = real(PS_avg_fwd./(1i+I_avg_fwd./V));

PS = PS-PS(:,end);

%% Plot data image
figure;
axes;

imagesc(Iset,V,PS);
% imagesc(Iset,V,PS./Iset');

colormap(jet(64));
set(gca, 'YDir', 'normal');
xlabel('Iset (nA)');
ylabel('E (meV)');
title('Iset spectroscopy');

yyaxis right;
plot(Iset,Z,'w','LineWidth',2);
ylabel('Z (nm)');
yyaxis left;

%% Plot a few curves
mid_idx = int16(round(length(datfiles)/2));
figure;
axes;
hold on;

plot(V,PS(:,1)./Iset(1),'LineWidth',2,...
    'DisplayName',['Iset = ' num2str(Iset(1)) ' nA']);
plot(V,PS(:,mid_idx)./Iset(mid_idx),'LineWidth',2,...
    'DisplayName',['Iset = ' num2str(Iset(mid_idx)) ' nA']);
plot(V,PS(:,end)./Iset(end),'LineWidth',2,...
    'DisplayName',['Iset = ' num2str(Iset(end)) ' nA']);

hold off;
box on;
xlabel('E (meV)');
ylabel('dI/dV/Iset (a.u.)');
title('Iset spectroscopy');
legend;

%% Plot data waterfall
N = size(PS,2);
cmap = jet(N);

figure;
hold on;
for i = 1:N-9
    plot(V,PS(:,i)+(i-1)/10,'LineWidth',2,'Color',cmap(i,:));
end
hold off;
box on;
xlabel('E (meV)');
ylabel('dI/dV (a.u.) + offset');
title('Iset spectroscopy');