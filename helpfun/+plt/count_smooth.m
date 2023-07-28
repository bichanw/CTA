function count_smooth(spk_count)
% plot smoothed spike count

% zscore first
spk_count = zscore(spk_count,[],2);

% currently don't know what binning parameters i should use
ax = np;
imagesc(ops.posterior_t,1:size(spk_count,1),spk_count);
xlim([600 700]);
ax.CLim = [-2 2];
ef;

% smooth by convolving with a half-gaussian
% calculate kernel
kernel_width = p.Results.kernel_width; % sigma
n_points = 3*kernel_width;
v = 1:n_points;
v = exp(-v.^2 / (2*kernel_width^2)) ./ (kernel_width * sqrt(2*pi));
v = v./v(1);

% assign spikes
spk_mat = zeros(n_cells,numel(t_mat));
for icell = 1:n_cells
	spk_mat(icell,ismember(t_mat,uint32(spikes{icell}*1e3))) = 1;
	spk_mat(icell,:) = conv(spk_mat(icell,:),v,'same');
end