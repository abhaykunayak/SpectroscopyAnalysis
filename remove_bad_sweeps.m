function [LS_clean, I_clean] = remove_bad_sweeps(DATA, sig)
% Removes bad sweeps

LS_sweeps = 5:3:size(DATA,4);
LS = DATA(:,:,:,5:3:end);
I = squeeze(DATA(:,:,:,4:3:end));

LS_avg = mean(LS, 4);
LS_std = std(LS, 0, 4);
LS_diff = DATA(:,:,:,5:3:end) - repmat(LS_avg, 1, 1, 1, numel(LS_sweeps));
LS_select = abs(LS_diff)<LS_std*sig;
LS_clean = mean(LS.*LS_select, 4);
I_clean = mean(I.*LS_select, 4);

end