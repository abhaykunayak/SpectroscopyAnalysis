function [Z_drift_corrected, LS_drift_corrected] = drift_correct(Z_subtracted_linear_fit, LS)

%% Plot data
figure;
axes;
imagesc(Z_subtracted_linear_fit);
% imagesc(Z_drift_corrected);
set(gca, 'YDir', 'Normal');
axis tight;
colormap(gray);
caxis auto

%% Gather points
alignment_points = round(ginput);
alignment_points(:,1) = alignment_points(:,1)-min(alignment_points(:,1))+1;

%% Shift data
Z_drift_corrected = zeros(size(Z_subtracted_linear_fit));
LS_drift_corrected = zeros(size(LS));

for i=1:size(alignment_points,1)
    line_no = alignment_points(i,2);
    Z_drift_corrected(line_no,:) = [Z_subtracted_linear_fit(line_no,alignment_points(i,1):end) zeros(1, alignment_points(i,1)-1)];
    LS_drift_corrected(line_no,:,:) = [LS(line_no,alignment_points(i,1):end,:) zeros(1, alignment_points(i,1)-1,size(LS,3))];
end

end