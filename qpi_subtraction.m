function [] = qpi_subtraction(ax1, ax2, himrect1, himrect2, en_shift_1, en_shift_2)

if nargin<5
    en_shift_1 = 0;
    en_shift_2 = 0;
end

% Plot subtraction figure
figure;
ax = axes;
im = imagesc(0,0,0);
xlabel('q_{x} (nm^{-1})');
ylabel('E (eV)');
title('Subtracted FFT');
axis tight;
set(gca,'YDir','normal');
colormap();
colorbar();
caxis(ax, [-1 1]);
addNewPositionCallback(himrect1, @(p) update_pos(p,himrect2,ax1,ax2,im,...
    en_shift_1,en_shift_2));
end

function [] = update_pos(p,himrect2,ax1,ax2,im,en_shift_1,en_shift_2)
set_xpos(p,himrect2);
update_plot(ax1, ax2, im, en_shift_1, en_shift_2);
end

function [] = set_xpos(p,himrect2)
p_old = getPosition(himrect2);
p_old(3) = p(3);
setPosition(himrect2,p_old)
end

function [] = update_plot(ax1, ax2, im, en_shift_1, en_shift_2)
% Energy shift correction
ax1.Children.CData = energy_offset(ax1.Children.YData, ...
    ax1.Children.CData, en_shift_1);
ax2.Children.CData = energy_offset(ax2.Children.YData, ...
    ax2.Children.CData, en_shift_2);

set(im, 'XData', ax1.Children.XData);
set(im, 'YData', ax1.Children.YData);
set(im, 'CData', ax1.Children.CData-1.5.*ax2.Children.CData);
end

function [ DATA_shift ] = energy_offset(V, DATA, V_I0)
%energy_offset shifts a 2D QPI/FFT DATA along the energy axis, given the
%bias vector and the bias at I=0.  

dV=abs(mean(diff(V)));
shift_ind=round(V_I0/dV);
DATA_shift=circshift(DATA,shift_ind,1);

end

