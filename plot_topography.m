function fig = plot_topography(X,Y,Z)
fig = figure;
imagesc(X, Y, Z);
xlabel('X (nm)','FontSize',14);
ylabel('Y (nm)','FontSize',14);
title('Topography','FontSize',16);
axis tight;
set(gca,'YDir','normal');
colormap('gray');
colorbar();
caxis auto;
end