function rate_smooth(spikes,toi,ax,varargin)

% parse input
p = inputParser;
addParameter(p,'kernel_width',25);
parse(p,varargin{:});


% intitiation
n_cells = numel(spikes);
t_mat = uint32((toi(1)*1e3):(toi(2)*1e3));
spk_mat = zeros(n_cells,numel(t_mat));


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

% bin for easier plotting
[R,~,~,tbin] = running_average(t_mat,spk_mat,10,10);


% plot
% ax = np;
imagesc(ax,tbin / 1000,1:n_cells,(R ./ max(R,[],1))');
% imagesc(ax,tbin / 1000,1:n_cells,R');
% imagesc(ax,tbin / 1000,1:n_cells,zscore(R,0,1)');
set(ax,'YDir','reverse','XLim',toi,'YLim',[0 n_cells]);
colorbar;
colormap(cbrewer2('Greys'));
set(gcf,'Position',[0 0 diff(toi)*8 n_cells*2]);
