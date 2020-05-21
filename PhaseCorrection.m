theta=2.2; 
PS_avg_X_rot = cosd(theta)*mean(PS_avg_X,2)-sind(theta)*mean(PS_avg_Y,2);
PS_avg_Y_rot = sind(theta)*mean(PS_avg_X,2)+cosd(theta)*mean(PS_avg_Y,2);

figure;
hold on;
plot(V, mean(PS_avg_X,2));
plot(V, mean(PS_avg_Y,2));

plot(V, PS_avg_X_rot);
plot(V, PS_avg_Y_rot);
hold off;