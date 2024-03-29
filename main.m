clear all; clc; close all; addpath('helpfun');

% ----- individual cell plotting ----- %
% % parameters
% plt_name = 'FR';
% extension = 'png';


sessions; % initialize sessions to run
for iSession = 1:numel(Sessions.subject)
% for iSession = 2


	% load data
	data = load_data.all(Sessions.session(iSession),Sessions.subject{iSession});
	
	% parameter section
	params;

	% ----- individual cell level plot ----- %
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


	% ----- session level plot ----- %
	% classifier.plt.psth(data,struct('tp',[0 1],'novel_vs_fam',struct('n_sig',15)));
	% classifier.plt.raster(data,struct('tp',[0 1],'novel_vs_fam',struct('n_sig',15)));
	% classifier.plt.cv(data,ops,sprintf('all_%s_',ops.mnr.penalty));
	% ops = struct('novel_vs_fam',struct('n_sig',15),'posterior_method',2);classifier.plt.posterior_raster(data,[],ops,'chris_thresholded_'); % plot posterior
	% ops.tp = [0 1]; classifier.plt.dprime(data,ops);
	% classifier.plt.avg_firing(data,ops);
	% classifier.plt.slow_firing(data,ops,'');
	

	% select by brain region
	% roi = {'CEA'};
	% ops.exclude_id = getOr(ops,'exclude_id')' | classifier.select_cells.by_brain_region(data,roi);

	% % for lambda = [1e-4 1e-2 1]
	for lambda = 1
		ops.mnr.lambda = lambda;
		prefix = sprintf('amp%d_%dcat_%s_%s_%.1e_',ops.amplitude_cutoff,numel(ops.amplitude_cat),ops.zscore_method,ops.mnr.penalty,ops.mnr.lambda);
		% prefix = 'testremove_';

		% classifier.plt.posterior_raster(data,[],ops,prefix); % plot posterior
		% [ops,to_save(iSession)] = classifier.plt.slow_firing(data,ops,prefix);
		ops.posterior_t_edges = data.video(1):0.15:data.video(end); [ops,resp(iSession,:,:,:),resp_err(iSession,:,:,:),tbin,raw_data(iSession,:)] = classifier.plt.laser_posterior(data,ops,prefix);
		% ops.posterior_t_edges = data.video(1):0.15:data.video(end); [ops,resp(iSession,:,:,:),resp_err(iSession,:,:,:),tbin,raw_data(iSession,:,:)] = classifier.plt.drinking_posterior(data,ops,prefix);
		% classifier.plt.examine_coef(data,ops,prefix);
		% classifier.plt.brain_region(data,ops,prefix); 
		% return
	end
	% return
	fprintf('%d session, front is %s\n',iSession,data.ports.front);
	% to_save(iSession).cv_results = classifier.plt.cv(data,ops,prefix);
	% return
end

% save processing setting
saveops(ops);
% return

if exist('to_save')
	% save(sprintf('mat/decoders_peaks_bin%d_step%d.mat',to_save(1).bin_width,to_save(1).step_size),'to_save');
	save(sprintf('mat/%s.mat',datestr(now,'YYmmdd')),'to_save');
end
if exist('resp')
	save(sprintf('mat/%s.mat',datestr(now,'YYmmdd')),'resp','resp_err','tbin','raw_data');
	% save('figures/0.15.mat','resp','resp_err','tbin','raw_data');
end

% append_script(['results/' prefix],true);
% append_script(['results/novel_vs_fam_' prefix],true);
% fprintf('Finish running %s\n',prefix);
return




% examine 
tmp = reshape(catcell(Mdl.DistributionParameters),2,3,110);
ax = np; imagesc(squeeze(tmp(1,:,:)));colorbar;ef;

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


