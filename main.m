clear all; clc; close all; addpath('helpfun');

% ----- individual cell plotting ----- %
% parameters
plt_name = 'FR';
extension = 'png';


sessions; % initialize sessions to run
for iSession = 1:numel(Sessions.subject)

	% load data
	data = load_data.all(Sessions.session(iSession),Sessions.subject{iSession});

	% parameter section
	params;

	% cell selection
	ops = classifier.select_cells.sig_resp(data,ops);
	% ops = struct('exclude_id',false(1,numel(data.spikes))); % all cells

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
	% ops = struct('tp',[0 1],'novel_vs_fam',struct('n_sig',15)); classifier.plt.psth(data,ops);
	% ops = struct('tp',[0 1],'novel_vs_fam',struct('n_sig',15)); classifier.plt.raster(data,struct());
	% ops.if_exclude_0var = true; ops.fake0 = 0;classifier.plt.cv(data,ops,'pythonmnr');
	% ops = struct('novel_vs_fam',struct('n_sig',15),'posterior_method',2);classifier.plt.posterior_raster(data,[],ops,'chris_thresholded_'); % plot posterior
	% ops.tp = [0 1]; classifier.plt.dprime(data,ops);
	% classifier.plt.avg_firing(data,ops);
	
	% figure 2c
	ops.novel_vs_fam.n_sig = 15;
	ops = classifier.select_cells.novel_vs_fam(data,ops);
	ops.exclude_id = ~ismember(1:numel(data.spikes),ops.novel_vs_fam.ordered_id(1:ops.novel_vs_fam.ordered_div(end-2)));
	ops.exclude_method = {'novel_vs_fam'};

	for lambda = logspace(-3,3,3)
		ops.mnr.lambda = lambda;
		prefix = sprintf('30cells_%s_%.1e_',ops.mnr.penalty,ops.mnr.lambda);
		ops = classifier.plt.posterior_raster(data,[],ops,prefix); % plot posterior
		classifier.plt.cv(data,ops,prefix);
		classifier.plt.examine_coef(data,ops,prefix);
	end

end

% save processing setting
saveops(ops);



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





