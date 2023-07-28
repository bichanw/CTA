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
% final_type = ones(size(final_type)); % plot all black to see raw raster
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

return





% check julia data firing rate
	load('mat/julia.mat');
	[~,count,edges,toi] = running_average(t,[],1,0.25); 
	ax = np; plot(toi,count);
	data = load_data.all(datetime(2023,2,14),'280');
	h = classifier.plt.scatter_event(data,ax,300);
	ef;

% read julia data
	f = 'mat/100s_30cell_280.jld2';
	h5disp(f)
	h5disp(f,'/_types/')
	tmp = h5read(f,'/single_stored_object');


	f = 'mat/280_230214_assignments.jld2';
	final_assignments = h5read(f,'/final_assignments');
	final_events      = h5read(f,'/final_events');
	

	% assign spikes to event type
	[~,Locb] = ismember(final_assignments,final_events.assignment_id);

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



% distribution parameter check

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
