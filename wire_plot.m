function [] = wire_plot(data, xdata, run_ar, spacing, scaling)
figure;
axes;
hold on;

for i=run_ar
    p = plot(xdata, smooth(scaling.*data(i,:)+i./spacing, 5));
    p.LineWidth = 1.5;
    p.Color = [0 0 0];
end
hold off;
axis tight;
box on;

xlabel('Energy');
ylabel('dI/dV (au)');
title('MDC')
end