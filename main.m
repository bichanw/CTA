clear all; clc; close all; addpath('helpfun');

% ----- individual cell plotting ----- %
% parameters
plt_name = 'FR';
extension = 'png';


sessions; % initialize sessions to run
for iSession = 1:numel(Sessions.subject)
% for iSession = 4


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
		prefix = sprintf('novel_fam_diff_all_sepzscore_%s_%.1e_',ops.mnr.penalty,ops.mnr.lambda);
		% prefix = sprintf('%s_%s_%.1e_',roi{1},ops.mnr.penalty,ops.mnr.lambda);
		ops = classifier.plt.slow_firing(data,ops,prefix);
		classifier.plt.posterior_raster(data,[],ops,prefix); % plot posterior
		ops = classifier.plt.laser_posterior(data,ops,prefix);
		% classifier.plt.examine_coef(data,ops,prefix);
		% classifier.plt.brain_region(data,ops,prefix); 
	end
	% classifier.plt.cv(data,ops,prefix);
end

% save processing setting
saveops(ops);



return


% read julia data
	f = 'mat/100s_30cell_280.jld2';
	h5disp(f)
	h5disp(f,'/_types/')
	tmp = h5read(f,'/single_stored_object');


	f = 'mat/280_230214_assignments.jld2';
	final_assignments = h5read(f,'/final_assignments');
	final_events      = h5read(f,'/final_events');
	

	% assign spikes to event type
	[~,Locb]   = ismember(final_assignments,final_events.assignment_id);

	tmp = find(Locb~=0);
	final_type = final_assignments;
	final_type(tmp) = final_events.seq_type(Locb(tmp));


	% load spikes
	load('mat/julia.mat');

	% plot spikes
	ax = np;

	Colors = [0 0 0;1 0 0;0 0 1];

	% plot spikes
	cell_oi  = unique(clu);
	event_oi = unique(final_type);
	for ineuron = 1:numel(cell_oi)
		for ievent = 1:numel(event_oi)
			h = my_scatter(t(clu==cell_oi(ineuron)&final_type==event_oi(ievent)),ineuron,ax,...
							3,'MarkerEdgeColor','none','MarkerFaceColor',Colors(ievent,:),'MarkerFaceAlpha',0.1);
		end
	end

	% plot other events
	% ?
	data = load_data.all(datetime(2023,2,14),'280');
	h = classifier.plt.scatter_event(data,ax,numel(cell_oi)+1);

	for t_start = 300:100:900
		set(ax,'XLim',[0 100] + t_start);
		export_fig(sprintf('results/ppseq_%d.png',t_start),'-r300');
	end



% check julia firing rate
	id = unique(clu);
	dt = max(t)-min(t);

	fr = arrayfun(@(ii) sum(clu==ii)/dt, id);


% save data for julia
	data = load_data.all(datetime(2023,2,14),'280');
	params;

	[~,cell_order] = classifier.select_cells.novel_vs_fam(data,ops);
	ops.exclude_id = ~ismember(1:numel(data.spikes),cell_order.ordered_id(1:(cell_order.ordered_div(3))));
	[count,ops,events_oi] = classifier.count_spk.events(data,ops);
	cellfun(@(x) sum(x(:)) / size(x,1) / diff(ops.tp), count)
	
	clu = [];
	t   = [];
	toi = [300 1000];
	cell_oi = find(~ops.exclude_id);
	for ii = 1:sum(~ops.exclude_id)
		tmp = data.spikes{cell_oi(ii)}(data.spikes{cell_oi(ii)}>toi(1)&data.spikes{cell_oi(ii)}<toi(2));
		clu = [clu; ones(numel(tmp),1)*ii];
		t   = [t; tmp];
	end

	% saving
	filename = sprintf('%s_%s',data.subject,datestr(data.session,'YYmmdd'));
	save('mat/julia.mat','clu','t','cell_oi','filename','toi');


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



% gamma distribution
	m = 50; v = 10;
	b = m / v;
	a = m * b;
	x = 30:70;
	y = gampdf(x,a,1/b);
	ax = np; plot(x,y);ef;

% scaled inverse-chi-square simulation
	v = 6;
	t2 = 0.15;
	x = 0:0.01:1;
	f = (t2*v/2)^(v/2)/gamma(v/2) * x.^(-v/2-1) .* exp(-t2*v./(2*x));
	ax = np; plot(x,f);ef;


% 2d dirichlet simulation
	x = 0:0.01:1;
	conc_param = 10.0;
	y = (x - x.^2).^(conc_param-1);
	y = y ./ sum(y);
	ax = np; plot(x,y);export_fig tmp.png -r300;
