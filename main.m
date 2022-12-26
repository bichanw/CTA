% ----- plotting ----- %
% parameters
plt_name = 'FR_amp';

for iSession = 1:numel(Sessions.subject)

	% load data
	data = load_data.all(Sessions.session(iSession),Sessions.subject{iSession});
	data.amp = load_data.spike_amp(data); % load amplitude

	% plotting initiation
	cells_2_plt = 1:numel(data.spikes); 
	N = numel(cells_2_plt);
	[ax,r,c] = np(N);

	% plot cell by cell
	for ii = 1:N
		eval(['plt.' plt_name '(data,ii,ax(ii));'])
	end

	% figure setting
	eval(['set_fig.' plt_name '(data,ax(sub2ind([c r],1,r)));'])

	% save
	export_fig(sprintf('results/%s_%s_%s.png',plt_name,datestr(data.session,'YYmmdd'),data.subject),'-m3');
end