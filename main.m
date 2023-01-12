clear all; clc; close all; addpath('helpfun');

% ----- individual cell plotting ----- %
% parameters
plt_name = 'FR';
extension = 'png';


sessions; % initialize sessions to run
for iSession = 1:numel(Sessions.subject)

	% load data
	data = load_data.all(Sessions.session(iSession),Sessions.subject{iSession});

	% cell selection
	% ops = struct('exclude_id',false(1,numel(data.spikes))); % all cells
	ops = classifier.select_cells.sig_resp(data,struct());
	% return

	% ----- individual cells ----- %
	% % plotting initiation
	% cells_2_plt = 1:numel(data.spikes); 
	% N = numel(cells_2_plt);
	% N_in_batch = 100;
	% for ibatch = 1:ceil(N/N_in_batch)
	% 	[ax,r,c] = np(N_in_batch);

	% 	% plot cell by cell
	% 	for ii = 1:min([N_in_batch N-(ibatch-1)*N_in_batch])
	% 		eval(['plt.' plt_name '(data,ii+(ibatch-1)*N_in_batch,ax(ii));'])
	% 	end

	% 	% figure setting
	% 	eval(['set_fig.' plt_name '(ax(sub2ind([c r],1,r)),data);'])
		
	% 	% save
	% 	export_fig(sprintf('results/%s_%s_%s_%d.%s',plt_name,datestr(data.session,'YYmmdd'),data.subject,ibatch,extension));

	% end
	% % combine plots
	% if strcmp(extension,'pdf') append_script(sprintf('results/%s_%s_%s',plt_name,datestr(data.session,'YYmmdd'),data.subject)); end


	% session level processing
	ops = classifier.plt.posterior(data,[],ops); % plot posterior
	% ops.tp = [0 1]; classifier.plt.dprime(data,ops);


	% % each cell as a pdf?
	% cells_2_plt = 1:numel(data.spikes); 
	% for ii = cells_2_plt
	% 	ax = np;
	% 	plt.FR_amp(data,ii,ax);
	% 	set_fig.FR_amp(data,ax);
	% 	return
	% end


end

% save processing setting
saveops(ops);



return


'/jukebox/witten/Chris/matlab/cz/neuropixels-cta/calca001-20221130.mat'


% examine if amplitude change
	% first 5 min vs last 5 min?
	t_2_cmp = 5 * 60; % 300 s

	% number of spikes available?
	N = numel(data.spikes);
	n_spk = [arrayfun(@(ii) sum(data.spikes{ii}<data.video(1)+t_2_cmp), 1:N)',...
			 arrayfun(@(ii) sum(data.spikes{ii}>=data.video(end)-t_2_cmp), 1:N)'];

	% send report
	fprintf('Session %s_%s\n',datestr(data.session,'YYmmdd'),data.subject);
	fprintf('Start video frame %.1f s, end video frame %.1f s, time of interest %d sec.\n',data.video(1),data.video(end),t_2_cmp);
	arrayfun(@(ii) fprintf('Cell %d available spikes: %d, %d\n',ii,n_spk(ii,1),n_spk(ii,2)), 1:numel(data.spikes));





