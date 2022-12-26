% 001, 002, 277
data = load_data.all(datetime(2022,11,30),'001');
data.amp = load_data.spike_amp(data); % load amplitude



% ----- plotting ----- %
% parameters
plt_name = 'FR_amp';
cells_2_plt = 1:numel(data.spikes); %1:9;


% initiation
N = numel(cells_2_plt);
[ax,r,c] = np(N);

% plot
for ii = 1:N
	yyaxis(ax(ii),'right'); plt.amp(data,ii,ax(ii));
	yyaxis(ax(ii),'left');  plt.FR(data,ii,ax(ii));
end

% figure setting
set_fig.FR(data,ax(sub2ind([c r],1,r)));
set_fig.amp(ax(sub2ind([c r],1,r)));

% save
export_fig(sprintf('results/%s_%s_%s.png',plt_name,datestr(data.session,'YYmmdd'),data.subject),'-m3');
